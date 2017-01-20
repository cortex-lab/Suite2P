function ops1 = reg2P(ops)
%%
if getOr(ops, 'doRegistration', 1)
    disp('running rigid registration');
else
    disp('skipping registration, but assembling binary file');
end

numPlanes = length(ops.planesToProcess);

chunk_align        = getOr(ops, {'chunk_align'}, 1);
nplanes            = getOr(ops, {'nplanes'}, 1);
nchannels          = getOr(ops, {'nchannels'}, 1);
ichannel           = getOr(ops, {'gchannel'}, 1);
rchannel           = getOr(ops, {'rchannel'}, 2);
red_align          = getOr(ops, {'AlignToRedChannel'}, 0);
BiDiPhase          = getOr(ops, {'BiDiPhase'}, 0);
LoadRegMean        = getOr(ops, {'LoadRegMean'}, 0);
ops.RegFileBinLocation = getOr(ops, {'RegFileBinLocation'}, []);
ops.splitFOV           = getOr(ops, {'splitFOV'}, [1 1]);


ops.smooth_time_space = getOr(ops, 'smooth_time_space', []);
fs = ops.fsroot;

%% find the mean frame after aligning a random subset

indx = 0;
try
   IMG = loadFramesBuff(fs{1}(1).name, 1, 1, 1); 
catch
    error('could not find any tif or tiff, check your path');
end
[Ly, Lx, ~, ~] = size(IMG);
ops.Ly = Ly;
ops.Lx = Lx;
% split into subsets (for high scanning resolution recordings)
[xFOVs, yFOVs] = get_xyFOVs(ops);

if ops.doRegistration
    IMG = GetRandFrames(fs, ops);    
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
    
    ops1 = cell(numPlanes, size(xFOVs,2));
    for i = 1:numPlanes
        for j = 1:size(xFOVs,2)
            ops1{i,j} = align_iterative(single(squeeze(IMG(yFOVs(:,j),xFOVs(:,j),...
                ops.planesToProcess(i),:))), ops);
        end
    end
    
    if ops.showTargetRegistration   
        PlotRegMean(ops1,ops);
        drawnow
    end
    clear IMG
else
    for i = 1:numPlanes
        ops1{i} = ops;
        ops1{i}.mimg = zeros(Ly, Lx);
     end
end

%% open files for registration
fid = cell(numPlanes, size(xFOVs,2));
for i = 1:numPlanes
    for j = 1:size(xFOVs,2)
        ops1{i,j}.RegFile = fullfile(ops.RegFileRoot, ...
            sprintf('%s_%s_%s_plane%d.bin', ops.mouse_name, ops.date, ops.CharSubDirs, i + (j-1)*numPlanes));
        regdir = fileparts(ops1{i,j}.RegFile);
        if ~exist(regdir, 'dir')
            mkdir(regdir);
        end
        
        % open bin file for writing
        fid{i,j}              = fopen(ops1{i,j}.RegFile, 'w');
        ops1{i,j}.DS          = [];
        ops1{i,j}.CorrFrame   = [];
        ops1{i,j}.mimg1       = zeros(ops1{i,j}.Ly, ops1{i,j}.Lx);
    end
end


%%
tic
% if two consecutive files have as many bytes, they have as many frames
nbytes = 0;
for k = 1:length(fs)
    for i = 1:numel(ops1)
         ops1{i}.Nframes(k)     = 0;
         ops1{i}.badframes(sum(ops1{i}.Nframes) + 1)   = true;
    end
    
    iplane0 = 1:1:ops.nplanes;
    for j = 1:length(fs{k})
        if abs(nbytes - fs{k}(j).bytes)>1e3
            nbytes = fs{k}(j).bytes;
            nFr = nFramesTiff(fs{k}(j).name);
        end
        
        iplane0 = mod(iplane0-1, numPlanes) + 1;
        if red_align
            ichanset = [rchannel; nFr; nchannels];
        else
            ichanset = [ichannel; nFr; nchannels];
        end
        data = loadFramesBuff(fs{k}(j).name, ichanset(1), ichanset(2), ...
            ichanset(3), ops.temp_tiff);
        
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
            [dsall, ops1] = GetRegOffsets(data, j, iplane0, ops, ops1, yFOVs, xFOVs);
            
            % if aligning by the red channel, data needs to be reloaded as the
            % green channel
            if red_align
%                 nFr = nFramesTiff(fs{k}(j).name);
                if mod(nFr, nchannels) ~= 0
                    fprintf('  WARNING: number of frames in tiff (%d) is NOT a multiple of number of channels!\n', j);
                end
                ichanset = [ichannel; nFr; nchannels];
                data = loadFramesBuff(ops.temp_tiff, ichanset(1), ichanset(2), ichanset(3), ops.temp_tiff);               
            end
            
            dreg = RegMovie(data, ops1, dsall, yFOVs, xFOVs);
        else
            dreg = data;
        end
        % write dreg to bin file+
        for i = 1:numPlanes
            ifr0 = iplane0(ops.planesToProcess(i));
            indframes = ifr0:nplanes:size(data,3);
            for l = 1:size(xFOVs,2)
                dwrite = dreg(yFOVs(:,l),xFOVs(:,l),indframes);
                fwrite(fid{i,l}, dwrite, class(data));
                
                ops1{i,l}.Nframes(k) = ops1{i,l}.Nframes(k) + size(dwrite,3);
                ops1{i,l}.mimg1 = ops1{i,l}.mimg1 + sum(dwrite,3);
            end
        end
        
        if rem(j,5)==1
            fprintf('Set %d, tiff %d done in time %2.2f \n', k, j, toc)            
        end
        
        iplane0 = iplane0 - nFr/nchannels;
    end    
end
for i = 1:numel(ops1)
    ops1{i}.mimg1 = ops1{i}.mimg1/sum(ops1{i}.Nframes);
    
    ops1{i}.badframes(1+sum(ops1{i}.Nframes)) = false;
    ops1{i}.badframes(1+sum(ops1{i}.Nframes)) = [];
end
%%
for i = 1:numPlanes    
    fclose(fid{i});
    fid{i}           = fopen(ops1{i}.RegFile, 'r');
    
    if ~isempty(ops.RegFileTiffLocation)
        ops1{i} = write_reg_to_tiff(fid{i}, ops1{i}, i);
    end    
    if ~isempty(ops.nimgbegend) && ops.nimgbegend>0
        ops1{i} = getBlockBegEnd(fid{i}, ops1{i}); % get mean of first and last frames in block (to check for drift)
    end
    if ~isempty(ops.RegFileBinLocation)
        str = sprintf('%d_',ops1{i}.expts);
        str(end) = [];
        folder = fullfile(ops1{i}.RegFileBinLocation, ops1{i}.mouse_name, ...
            ops1{i}.date, str);
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
for i = 1:numel(ops1)
    if ops.doRegistration
        % determine bad frames
        badi                    = getOutliers(ops1{i}); 
        ops1{i}.badframes(badi) = true;
        
        minDs = min(ops1{i}.DS(~ops1{i}.badframes, [1 2]), [], 1);
        maxDs = max(ops1{i}.DS(~ops1{i}.badframes, [1 2]), [], 1);
        disp([minDs(1) maxDs(1) minDs(2) maxDs(2)])
        if BiDiPhase>0
            maxDs(2) = max(1+BiDiPhase, maxDs(2));
        elseif BiDiPhase<0
            minDs(2) = min(BiDiPhase, minDs(2));
        end
        
        ops1{i}.yrange = ceil(1 + maxDs(1)) : floor(ops1{i}.Ly+minDs(1));
        ops1{i}.xrange = ceil(1 + maxDs(2)) : floor(ops1{i}.Lx+minDs(2));
    else
        ops1{i}.yrange = 1:Ly;
        ops1{i}.xrange = 1:Lx;
    end
    
    savepath = sprintf('%s/', ops.ResultsSavePath);
    
    if ~exist(savepath, 'dir')
        mkdir(savepath)
    end   
end

save(sprintf('%s/regops_%s_%s.mat', ops.ResultsSavePath, ...
    ops.mouse_name, ops.date),  'ops1')


%save(sprintf('%s/F_%s_%s_plane%d.mat', ops.ResultsSavePath, ...
 %   ops.mouse_name, ops.date, iplane), 'ops')

%%
