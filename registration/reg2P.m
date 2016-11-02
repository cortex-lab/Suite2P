function ops1 = reg2P(ops)
%%
numPlanes = length(ops.planesToProcess);

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


ops.smooth_time_space = getOr(ops, 'smooth_time_space', []);
fs = ops.fsroot;

%% find the mean frame after aligning a random subset
ntifs = sum(cellfun(@(x) numel(x), fs));
nfmax = ceil(ops.NimgFirstRegistration/ntifs);
if nfmax>=2000
    nfmax = 1999;
end
nbytes = fs{1}(1).bytes;
nFr = nFrames(fs{1}(1).name);
Info0 = imfinfo(fs{1}(1).name);
Ly = Info0(1).Height;
Lx = Info0(1).Width;
ops.Lx = Lx;
ops.Ly = Ly;

[xFOVs, yFOVs] = get_xyFOVs(ops); 

indx = 0;
IMG = zeros(Ly, Lx, nplanes, ops.NimgFirstRegistration, 'int16');

if ops.doRegistration
    for k = 1:length(ops.SubDirs)
        iplane0 = 1;
        for j = 1:length(fs{k})
            if abs(nbytes - fs{k}(j).bytes)>1e3
                nbytes = fs{k}(j).bytes;
                nFr = nFrames(fs{k}(j).name);
            end
            if nFr<(nchannels*nplanes*nfmax + nchannels*nplanes)
                continue;
            end
            
            iplane0 = mod(iplane0-1, nplanes) + 1;
            offset = 0;
            if j==1
                offset = nchannels*nplanes;
            end
            if red_align
                ichanset = [offset + nchannels*(iplane0-1) + [rchannel;...
                    nchannels*nplanes*nfmax]; nchannels];
            else
                ichanset = [offset + nchannels*(iplane0-1) + [ichannel;...
                    nchannels*nplanes*nfmax]; nchannels];
            end
            iplane0 = iplane0 - nFr/nchannels;
            data = loadFramesBuff(fs{k}(j).name, ichanset(1),ichanset(2), ichanset(3));
            data = reshape(data, Ly, Lx, nplanes, []);
            
            if BiDiPhase
                yrange = 2:2:Ly;
                if BiDiPhase>0
                    data(yrange, (1+BiDiPhase):Lx,:,:) = data(yrange, 1:(Lx-BiDiPhase),:,:);
                else
                    data(yrange, 1:Lx+BiDiPhase,:,:) = data(yrange, 1-BiDiPhase:Lx,:,:);
                end
            end
            IMG(:,:,:,indx+(1:size(data,4))) = data;
            indx = indx + size(data,4);
            
        end
    end
    IMG =  IMG(:,:,:,1:indx);
    
    ops1 = cell(numPlanes, size(xFOVs,2));
    for i = 1:numPlanes
        for j = 1:size(xFOVs,2)
            ops1{i,j} = align_iterative(single(squeeze(IMG(yFOVs(:,j),xFOVs(:,j),...
                ops.planesToProcess(i),:))), ops);
        end
    end
    
    if ops.showTargetRegistration
        nRows = floor(sqrt(numel(ops1)));
        nColumns = ceil(numel(ops1)/nRows);
        figure('Name', 'Registration Target Frames', ...
            'Position', [50 50 nColumns*500 nRows*500])
        i0 = 0;
        for i = 1:numPlanes
            for j = 1:size(xFOVs,2)
                i0 = i0+1;
                subplot(nRows, nColumns, i0);
                imagesc(ops1{i,j}.mimg);
                colormap('gray')
                axis equal tight
                title(sprintf('Plane %d, %s_%s', i, ops.mouse_name, ops.date))
            end
        end
        
        drawnow
    end
    
    clear IMG
else
    for i = 1:numPlanes
        ops1{i} = ops;
        ops1{i}.mimg = zeros(Ly, Lx);
        ops1{i}.Ly   = Ly;
        ops1{i}.Lx   = Lx;
     end
end
%%
fid = cell(numPlanes, size(xFOVs,2));
for i = 1:numPlanes
    for j = 1:size(xFOVs,2)
        ops1{i,j}.RegFile = fullfile(ops.RegFileRoot, sprintf('tempreg_plane%d.bin', i + (j-1)*numPlanes));
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
for k = 1:length(fs)
    for i = 1:numel(ops1)
         ops1{i}.Nframes(k)  = 0;
    end
    
    iplane0 = 1:1:ops.nplanes;
    for j = 1:length(fs{k})
        iplane0 = mod(iplane0-1, numPlanes) + 1;
       
        if red_align
            %                 ichanset = [nchannels*(iplane0-1)+rchannel; nFr; nchannels];
            ichanset = [rchannel; NaN; nchannels];
        else
            % using NaN instead of nFrames, we want to ask for nFrames
            % after copying the file locally (way faster)
            ichanset = [ichannel; NaN; nchannels];
        end
        data = loadFramesBuff(fs{k}(j).name, ichanset(1), ichanset(2), ichanset(3), ops.temp_tiff);
        nFr = size(data, 3)*nchannels;
        
        if BiDiPhase
            yrange = 2:2:Ly;
            if BiDiPhase>0
                data(yrange, (1+BiDiPhase):Lx,:) = data(yrange, 1:(Lx-BiDiPhase),:);
            else
                data(yrange, 1:Lx+BiDiPhase,:) = data(yrange, 1-BiDiPhase:Lx,:);
            end
        end
        
        if ops.doRegistration
            % get the registration offsets
            dsall = zeros(size(data,3), 2, size(xFOVs,2));
            for i = 1:numPlanes
                ifr0 = iplane0(ops.planesToProcess(i));
                indframes = ifr0:nplanes:size(data,3);

                for l = 1:size(xFOVs,2)
                    dat = data(yFOVs(:,l),xFOVs(:,l),indframes);
                    if ~isempty(ops.smooth_time_space)
                        dat = smooth_movie(dat, ops); 
                    end
                    [ds, Corr]  = registration_offsets(dat, ops1{i,l}, 0);
                    dsall(indframes,:, l)  = ds;
                    % collect ds
                    if j==1
                        ds(1,:,:) = 0;
                    end
                    ops1{i,l}.DS          = cat(1, ops1{i,l}.DS, ds);
                    ops1{i,l}.CorrFrame   = cat(1, ops1{i,l}.CorrFrame, Corr);
                end
            end
            
            % if aligning by the red channel, data needs to be reloaded as the
            % green channel
            if red_align
                nFr = nFramesTiff(fs{k}(j).name);
                if mod(nFr, nchannels) ~= 0
                    fprintf('  WARNING: number of frames in tiff (%d) is NOT a multiple of number of channels!\n', j);
                end
                ichanset = [ichannel; nFr; nchannels];
%                 data = img.loadFrames(fs{k}(j).name, ichanset(1), ichanset(2), ichanset(3));           
                % providing the last argument ops.temp_tiff to let the
                % function know that this is already a local file (will load faster)
                data = loadFramesBuff(ops.temp_tiff, ichanset(1), ichanset(2), ichanset(3), ops.temp_tiff);               
            end
            
            ix0 = 0;
            Nbatch = 1000;
            dreg = zeros(size(data), 'like', data);            
            while ix0<size(data,3)
                indxr = ix0 + (1:Nbatch);
                indxr(indxr>size(data,3)) = [];
                for l = 1:size(xFOVs,2)
                    dreg(yFOVs(:,l), xFOVs(:,l), indxr)        = ...
                        register_movie(data(yFOVs(:,l), xFOVs(:,l), indxr), ops1{1,l}, dsall(indxr,:,l));
                end
                ix0 = ix0 + Nbatch;
            end
            
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
        minDs = min(ops1{i}.DS(2:end, [1 2]), [], 1);
        maxDs = max(ops1{i}.DS(2:end, [1 2]), [], 1);
        disp([minDs(1) maxDs(1) minDs(2) maxDs(2)])
        if BiDiPhase>0
            maxDs(2) = max(1+BiDiPhase, maxDs(2));
        elseif BiDiPhase<0
            minDs(2) = min(BiDiPhase, minDs(2));
        end
        
        ops1{i}.yrange = ceil(maxDs(1)):floor(ops1{i}.Ly+minDs(1));
        ops1{i}.xrange = ceil(maxDs(2)):floor(ops1{i}.Lx+minDs(2));
    else
        ops1{i}.yrange = 1:Ly;
        ops1{i}.xrange = 1:Lx;
    end
    
    savepath = sprintf('%s/', ops.ResultsSavePath);
    
    if ~exist(savepath, 'dir')
        mkdir(savepath)
    end
    ops = ops1{i};
    save(sprintf('%s/regops_%s_%s_plane%d.mat', ops.ResultsSavePath, ...
        ops.mouse_name, ops.date, i),  'ops')
end


%save(sprintf('%s/F_%s_%s_plane%d.mat', ops.ResultsSavePath, ...
 %   ops.mouse_name, ops.date, iplane), 'ops')

%%
