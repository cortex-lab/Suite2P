% registers frames using block registration in Y
function ops1 = blockReg2P(ops)
%%
if ops.doRegistration
    disp('running non-rigid registration');
else
    disp('skipping registration, but assembling binary file');
end

% default is 8 blocks in the Y direction (1/6 pixels each)
ops.numBlocks      = getOr(ops, {'numBlocks'}, [8 1]);
numBlocks          = ops.numBlocks;
numPlanes = length(ops.planesToProcess);
ops.numPlanes = numPlanes;
chunk_align        = getOr(ops, {'chunk_align'}, 1);
nplanes            = getOr(ops, {'nplanes'}, 1);
nchannels          = getOr(ops, {'nchannels'}, 1);
ichannel           = getOr(ops, {'gchannel'}, 1);
rchannel           = getOr(ops, {'rchannel'}, 2);
red_align   = getOr(ops, {'AlignToRedChannel'}, 0);
BiDiPhase          = getOr(ops, {'BiDiPhase'}, 0);
LoadRegMean        = getOr(ops, {'LoadRegMean'}, 0);
ops.RegFileBinLocation = getOr(ops, {'RegFileBinLocation'}, []);
ops.splitFOV           = getOr(ops, {'splitFOV'}, [1 1]);
ops.RegFileTiffLocation = getOr(ops, {'RegFileTiffLocation'}, []);
ops.dobidi         = getOr(ops, {'dobidi'}, 1);

fs = ops.fsroot;

try
   IMG = loadFramesBuff(fs{1}(1).name, 1, 1, 1); 
catch
    error('could not find any tif or tiff, check your path');
end

% find the mean frame after aligning a random subset
if ops.doRegistration
    [IMG] = GetRandFrames(fs, ops);    
    [Ly, Lx, ~, ~] = size(IMG);
    ops.Ly = Ly;
    ops.Lx = Lx;
    ops = MakeBlocks(ops);
    fprintf('--- using %d blocks in Y\n', numBlocks(1));
    fprintf('--- %d pixels/block; avg pixel overlap = %d pixels\n', ...
        round(ops.blockFrac(1)*Ly), ops.pixoverlap(1) );

    if numBlocks(2) > 1
        fprintf('--- using %d blocks in X\n', numBlocks(2));
        fprintf('--- %d pixels/block; avg pixel overlap = %d pixels\n', ...
            round(ops.blockFrac(2)*Lx),  ops.pixoverlap(2));
    end

    
    % compute phase shifts from bidirectional scanning
    if ops.dobidi
        BiDiPhase = BiDiPhaseOffsets(IMG);
    else
        BiDiPhase = 0;
    end
    fprintf('bi-directional scanning offset = %d pixels\n', BiDiPhase);
   
    if abs(BiDiPhase) > 0
        yrange = 2:2:Ly;
        if BiDiPhase>0
            IMG(yrange,(1+BiDiPhase):Lx,:,:) = IMG(yrange, 1:(Lx-BiDiPhase),:,:);
        else
            IMG(yrange,1:Lx+BiDiPhase,:,:)   = IMG(yrange, 1-BiDiPhase:Lx,:,:);
        end
    end
    for i = 1:numPlanes
        ops1{i} = blockAlignIterative(squeeze(IMG(:,:,ops.planesToProcess(i),:)), ops);
    end
   if ops.showTargetRegistration
       PlotRegMean(ops1,ops);
   end
else
    [Ly, Lx] = size(IMG);
    for i = 1:numPlanes
        ops1{i} = ops;
        ops1{i}.mimg = zeros(Ly, Lx);
        ops1{i}.Ly   = Ly;
        ops1{i}.Lx   = Lx;
    end
end
clear IMG


% prepare individual options files and open binaries
for i = 1:numPlanes
    ops1{i}.RegFile = fullfile(ops.RegFileRoot, sprintf('tempreg_plane%d.bin',...
        ops.planesToProcess(i)));
    regdir = fileparts(ops1{i}.RegFile);
    if ~exist(regdir, 'dir')
        mkdir(regdir);
    end    
    % open bin file for writing
    fid{i}             = fopen(ops1{i}.RegFile, 'w');
    ops1{i}.DS          = [];
    ops1{i}.CorrFrame   = [];
    ops1{i}.mimg1       = ops1{i}.mimg;
    for ib = 1:numBlocks(1)*numBlocks(2)
        ops1{i}.mimgB{ib} = ops1{i}.mimg(ops1{i}.yBL{ib}, ops1{i}.xBL{ib});
    end
end
nbytes = fs{1}(1).bytes;
nFr = nFramesTiff(fs{1}(1).name);

% compute registration offsets from mean img for each frame
xyValid = true(Ly, Lx);
tic
for k = 1:length(fs)
    for i = 1:numPlanes
         ops1{i}.Nframes(k)  = 0;
    end
    iplane0 = 1:1:ops.nplanes;
    for j = 1:length(fs{k})
        iplane0 = mod(iplane0-1, numPlanes) + 1;
        nFr = nFramesTiff(fs{k}(j).name);
        if red_align
            ichanset = [rchannel; nFr; nchannels];
        else
            ichanset = [ichannel; nFr; nchannels];
        end
        data = loadFramesBuff(fs{k}(j).name, ichanset(1), ichanset(2), ichanset(3), ops.temp_tiff);

        if abs(BiDiPhase) > 0
            yrange = 2:2:Ly;
            if BiDiPhase>0
                data(yrange, (1+BiDiPhase):Lx,:) = data(yrange, 1:(Lx-BiDiPhase),:);
            else
                data(yrange, 1:Lx+BiDiPhase,:) = data(yrange, 1-BiDiPhase:Lx,:);
            end
        end
        
        if ops.doRegistration
            % get the registration offsets
            [dsall, ops1] = GetBlockOffsets(data, j, iplane0, ops, ops1);
                    
            % if aligning by the red channel, data needs to be reloaded as the
            % green channel
            if red_align
                if mod(nFr, nchannels) ~= 0
                    fprintf('  WARNING: number of frames in tiff (%d) is NOT a multiple of number of channels!\n', j);
                end
                ichanset = [ichannel; nFr; nchannels];
%                 data = img.loadFrames(fs{k}(j).name, ichanset(1), ichanset(2), ichanset(3));                
                data = loadFramesBuff(ops.temp_tiff, ichanset(1), ichanset(2), ichanset(3));               
            end
            
            % register frames
            [dreg, xyValid] = BlockRegMovie(data, ops, dsall, xyValid);
       
        else
            dreg = data;
        end
        % write dreg to bin file+
        for i = 1:numPlanes
            ifr0 = iplane0(ops.planesToProcess(i));
            indframes = ifr0:nplanes:size(data,3);
            dwrite = dreg(:,:,indframes);
            fwrite(fid{i}, dwrite, class(data));
            
            ops1{i}.Nframes(k) = ops1{i}.Nframes(k) + size(dwrite,3);
            ops1{i}.mimg1 = ops1{i}.mimg1 + sum(dwrite,3);
        end
        
        if rem(j,5)==1
            fprintf('Set %d, tiff %d done in time %2.2f \n', k, j, toc)            
        end
        
        iplane0 = iplane0 - nFr/nchannels;
    end
    for i = 1:numPlanes
        ops1{i}.mimg1 = ops1{i}.mimg1/ops1{i}.Nframes(k);
    end    
end

%% write binary file to tiffs if ~isempty(ops.RegFileTiffLocation)
for i = 1:numPlanes    
    fclose(fid{i});
    fid{i}           = fopen(ops1{i}.RegFile, 'r');
    
    if ~isempty(ops.RegFileTiffLocation)
        ops1{i} = write_reg_to_tiff(fid{i}, ops1{i}, i);
        frewind(fid{i});
    end    
    if ~isempty(ops.RegFileBinLocation)
        folder = fullfile(ops1{i}.RegFileBinLocation, ops1{i}.mouse_name, ...
            ops1{i}.date);
        if ~exist(folder, 'dir')
            mkdir(folder)
        end
        fidCopy = fopen(fullfile(folder, sprintf('plane%d.bin', i)), 'w');
        sz = ops1{i}.Lx * ops1{i}.Ly;
        parts = ceil(sum(ops1{i}.Nframes) / 2000);
        for p = 1:parts
            toRead = 2000;
            if p == parts
                toRead = sum(ops1{i}.Nframes) - 2000 * (parts-1);
            end
            data = fread(fid{i},  sz*toRead, '*int16');
            fwrite(fidCopy, data, class(data));
        end
        fclose(fidCopy);
    end
    fclose(fid{i});
end
%%
% compute xrange, yrange
for i = 1:numPlanes
    if ops.doRegistration
%         badi                    = getOutliers(ops1{i}); 
%         ops1{i}.badframes(badi) = true;
        ops1{i}.badframes = false(sum(ops1{i}.Nframes), 1);
        
        %         minDs = min(min(ops1{i}.DS(2:end, :,:), [], 3), [], 1);
        %         maxDs = max(max(ops1{i}.DS(2:end, :,:), [], 3), [], 1);
        ops1{i}.yrange  = find(xyValid(:, ceil(Lx/2)));
        ops1{i}.xrange  = find(xyValid(ceil(Ly/2), :));
        
    else
        ops1{i}.yrange = 1:Ly;
        ops1{i}.xrange = 1:Lx;
    end
    
    savepath = sprintf('%s/', ops.ResultsSavePath);
    
    if ~exist(savepath, 'dir')
        mkdir(savepath)
    end
    
    mimg = zeros(size(ops1{i}.mimg1));
    for ib = 1:numBlocks
        mimg(ops1{i}.yBL{ib}, ops1{i}.xBL{ib}) = ops1{i}.mimgB{ib};
    end
    ops1{i}.mimg = mimg;
    
end

save(sprintf('%s/regops_%s_%s.mat', ops.ResultsSavePath, ...
        ops.mouse_name, ops.date),  'ops1')

%save(sprintf('%s/F_%s_%s_plane%d.mat', ops.ResultsSavePath, ...
 %   ops.mouse_name, ops.date, iplane), 'ops')
