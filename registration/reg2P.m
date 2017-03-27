function ops1 = reg2P(ops)
%%
if getOr(ops, 'doRegistration', 1)
    disp('running rigid registration');
else
    disp('skipping registration, but assembling binary file');
end

numPlanes = length(ops.planesToProcess);

nplanes            = getOr(ops, {'nplanes'}, 1);
nchannels          = getOr(ops, {'nchannels'}, 1);
ichannel           = getOr(ops, {'gchannel'}, 1);
rchannel           = getOr(ops, {'rchannel'}, 2);
red_align          = getOr(ops, {'AlignToRedChannel'}, 0);

ops.RegFileBinLocation = getOr(ops, {'RegFileBinLocation'}, []);
RegFileBinLocation = getOr(ops, {'RegFileBinLocation'}, []);
targetImage        = getOr(ops, {'targetImage'}, []); % specify experiment to generate target image from (useful if drift) 
alignTargetImages  = getOr(ops, {'alignTargetImages'}, false); % if true, align target images to each other
interpolateAcrossPlanes = getOr(ops, {'interpolateAcrossPlanes'}, false); %if true, similar looking planes will be averaged together to generate final movie
planesToInterpolate = getOr(ops, {'planesToInterpolate'}, 1:nplanes); % these planes will be considered for interpolation
alignAcrossPlanes  = getOr(ops, {'alignAcrossPlanes'}, false); % at each time point, frame will be aligned to best matching target image (from different planes)
ops.interpolateAcrossPlanes = getOr(ops, 'interpolateAcrossPlanes', 0);

ops.splitFOV           = getOr(ops, {'splitFOV'}, [1 1]);
ops.smooth_time_space  = getOr(ops, 'smooth_time_space', []);
ops.dobidi         = getOr(ops, {'dobidi'}, 1);
LoadRegMean        = getOr(ops, {'LoadRegMean'}, 0);


fs = ops.fsroot;

%% find the mean frame after aligning a random subset

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
    if ~getOr(ops, {'BiDiPhase'}, 0) && ops.dobidi
        BiDiPhase = BiDiPhaseOffsets(IMG);
    else
        BiDiPhase = getOr(ops, {'BiDiPhase'}, 0);
    end
    ops.BiDiPhase = BiDiPhase;
    fprintf('bi-directional scanning offset = %d pixels\n', BiDiPhase);
   
    if abs(BiDiPhase) > 0
        yrange = 2:2:Ly;
        if BiDiPhase>0
            IMG(yrange,(1+BiDiPhase):Lx,:,:) = IMG(yrange, 1:(Lx-BiDiPhase),:,:);
        else
            IMG(yrange,1:Lx+BiDiPhase,:,:)   = IMG(yrange, 1-BiDiPhase:Lx,:,:);
        end
    end
    
    % for each plane: align chosen frames to average to generate target image
    ops1 = cell(numPlanes, size(xFOVs,2));
    for i = 1:numPlanes
        for j = 1:size(xFOVs,2)
            ops1{i,j} = align_iterative(single(squeeze(IMG(yFOVs(:,j),xFOVs(:,j),...
                ops.planesToProcess(i),:))), ops);
        end
    end
    
    if alignTargetImages % align target images of all planes to each other
        % (reasonable only if interplane distance during imaging was small)
        ds = zeros(length(planesToInterpolate), 2, size(xFOVs,2));
        % planesToInterpolate should be indices of non-flyback planes
        % (only those should be aligned)
        for i = 2:length(planesToInterpolate) % align consecutive planes
            pl = planesToInterpolate(i);
            for j = 1:size(xFOVs,2)
                ds(i,:,j) = registration_offsets(ops1{pl,j}.mimg, ops1{pl-1,j}, 0);
            end
        end
        ds = cumsum(ds,1);
        ds = bsxfun(@minus, ds, mean(ds,1));
        images = zeros(Ly, Lx, length(planesToInterpolate), size(xFOVs,2));
        for i = 1:length(planesToInterpolate)
            for j = 1:size(xFOVs,2)
                images(:,:,i,j) = ops1{planesToInterpolate(i),j}.mimg;
            end
        end
        ds = reshape(permute(ds, [1 3 2]), [], 2);
        images = reshape(images, Ly, Lx, []);
        newTargets = register_movie(images, ops, ds);
        newTargets = reshape(newTargets, Ly, Lx, ...
            length(planesToInterpolate), size(xFOVs,2));
        for i = 1:length(planesToInterpolate)
            for j = 1:size(xFOVs,2)
                ops1{planesToInterpolate(i),j}.mimg = newTargets(:,:,i,j);
            end
        end
    end
    
    if ops.showTargetRegistration   
        PlotRegMean(ops1,ops);
        drawnow
    end
    clear IMG
else % no registration
    ops1 = cell(numPlanes, 1);
    for i = 1:numPlanes
        ops1{i} = ops;
        ops1{i}.mimg = zeros(Ly, Lx);
    end
end

%% open files for registration
fid = cell(numPlanes, size(xFOVs,2));
fidIntpol = [];
if ops.interpolateAcrossPlanes == 1 && ~isempty(RegFileBinLocation)
    fidIntpol = cell(numPlanes, size(xFOVs,2));
end
for i = 1:numPlanes
    for j = 1:size(xFOVs,2)
        ops1{i,j}.RegFile = fullfile(ops.RegFileRoot, ...
            sprintf('%s_%s_%s_plane%d.bin', ops.mouse_name, ops.date, ...
            ops.CharSubDirs, i + (j-1)*numPlanes));
        regdir = fileparts(ops1{i,j}.RegFile);
        if ~exist(regdir, 'dir')
            mkdir(regdir);
        end
        
        % open bin file for writing
        fid{i,j}              = fopen(ops1{i,j}.RegFile, 'w');
        ops1{i,j}.DS          = [];
        ops1{i,j}.CorrFrame   = [];
        ops1{i,j}.mimg1       = zeros(ops1{i,j}.Ly, ops1{i,j}.Lx);
        
        if interpolateAcrossPlanes && ~isempty(RegFileBinLocation)
            % open separate files for result after averaging across
            % neighbouring planes
            str = sprintf('%d_',ops1{i,j}.expts);
            str(end) = [];
            folder = fullfile(ops1{i,j}.RegFileBinLocation, ops1{i,j}.mouse_name, ...
                ops1{i,j}.date, str, 'interpolated');
            if ~exist(folder, 'dir')
                mkdir(folder);
            end
            fidIntpol{i,j} = fopen(fullfile(folder, ...
                sprintf('%s_%s_%s_plane%d.bin', ops.mouse_name, ops.date, ...
                ops.CharSubDirs, i + (j-1)*numPlanes)));
        end
    end
end

%%
tic
% if two consecutive files have as many bytes, they have as many frames
nbytes = 0;
for k = 1:length(fs)
    % initialize frame count
    for i = 1:numel(ops1)
         ops1{i}.Nframes(k)     = 0;
    end
    
    iplane0 = 1:1:ops.nplanes; % identity of planes for first frames in tiff file
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
        % only load frames of registration channel
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
            [dsall, ops1] = GetRegOffsets(data, k, j, iplane0, ops, ops1);
            
            if ~alignAcrossPlanes
                if red_align
                    % if aligning by the red channel, data needs to be reloaded as the
                    % green channel
                    %                 nFr = nFramesTiff(fs{k}(j).name);
                    if mod(nFr, nchannels) ~= 0
                        fprintf('  WARNING: number of frames in tiff (%d) is NOT a multiple of number of channels!\n', j);
                    end
                    ichanset = [ichannel; nFr; nchannels];
                    data = loadFramesBuff(ops.temp_tiff, ichanset(1), ichanset(2), ichanset(3), ops.temp_tiff);
                    if abs(BiDiPhase) > 0
                        yrange = 2:2:Ly;
                        if BiDiPhase>0
                            data(yrange, (1+BiDiPhase):Lx,:) = data(yrange, 1:(Lx-BiDiPhase),:);
                        else
                            data(yrange, 1:Lx+BiDiPhase,:) = data(yrange, 1-BiDiPhase:Lx,:);
                        end
                    end
                end
                % align the frames according to the registration offsets
                dreg = RegMovie(data, ops1, dsall, yFOVs, xFOVs);
            end
        else
            dreg = data;
        end
        
        if ~alignAcrossPlanes
            % write dreg to bin file+
            for i = 1:numPlanes
                ifr0 = iplane0(ops.planesToProcess(i));
                indframes = ifr0:nplanes:size(data,3);
                for l = 1:size(xFOVs,2)
                    dwrite = dreg(yFOVs(:,l),xFOVs(:,l),indframes);
                    fwrite(fid{i,l}, dwrite, class(data));
                    
                    ops1{i,l}.mimg1 = ops1{i,l}.mimg1 + sum(dwrite,3);
                end
            end
        end
        
        if rem(j,5)==1
            fprintf('Set %d, tiff %d done in time %2.2f \n', k, j, toc)            
        end
        
        iplane0 = iplane0 - nFr/nchannels;
    end    
end

if alignAcrossPlanes && ops.doRegistration % align each frame with the best matching target image
    [ops1, planesToInterp] = registerAcrossPlanes(ops1, ops, fid, fidIntpol);
end

for i = 1:numel(ops1)
    ops1{i}.mimg1 = ops1{i}.mimg1/sum(ops1{i}.Nframes);
    
    ops1{i}.badframes = false(1, size(ops1{i}.DS,1));
end
%%
for i = 1:numel(ops1)    
    fclose(fid{i});
    fid{i}           = fopen(ops1{i}.RegFile, 'r');
    
    if ~isempty(ops.RegFileTiffLocation)
        ops1{i} = write_reg_to_tiff(fid{i}, ops1{i}, i);
        frewind(fid{i});
    end    
    if ~isempty(ops.nimgbegend) && ops.nimgbegend>0
        ops1{i} = getBlockBegEnd(fid{i}, ops1{i}); % get mean of first and last frames in block (to check for drift)
        frewind(fid{i});
    end
    if ~isempty(ops.RegFileBinLocation)
        folder = fullfile(ops1{i}.RegFileBinLocation, ops1{i}.mouse_name, ...
            ops1{i}.date, ops.CharSubDirs);
        if ~exist(folder, 'dir')
            mkdir(folder)
        end
        fidCopy = fopen(fullfile(folder, ...
            sprintf('%s_%s_%s_plane%d.bin', ops.mouse_name, ops.date, ...
            ops.CharSubDirs, i)), 'w');
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
        fclose(fid{i});
        
        if ops.interpolateAcrossPlanes == 1 && ~isempty(RegFileBinLocation)
            fclose(fidIntpol{i});
            folder = fullfile(ops1{i}.RegFileBinLocation, ops1{i}.mouse_name, ...
                ops1{i}.date, ops.CharSubDirs, 'interpolated');
            filename = fullfile(folder, ...
                sprintf('%s_%s_%s_plane%d.bin', ops.mouse_name, ops.date, ...
                ops.CharSubDirs, i));
            if ismember(i, planesToInterp)
                fidCopy = fopen(ops1{i}.RegFile, 'w');
                fidOrig = fopen(filename, 'r');
                sz = ops1{i}.Lx * ops1{i}.Ly;
                parts = ceil(sum(ops1{i}.Nframes) / 2000);
                for p = 1:parts
                    toRead = 2000;
                    if p == parts
                        toRead = sum(ops1{i}.Nframes) - 2000 * (parts-1);
                    end
                    data = fread(fidOrig,  sz*toRead, '*int16');
                    fwrite(fidCopy, data, class(data));
                end
                fclose(fidCopy);
                fclose(fidOrig);
            else
                delete(filename)
            end
        end
    end
end
%%
% compute xrange, yrange
for i = 1:numel(ops1)
    if ops.doRegistration
        if size(ops1{i}.DS,3) > 1
            % determine bad frames
            ops = ops1{i};
            ops.CorrFrame(~ops.usedPlanes) = NaN;
            ops.CorrFrame = nanmean(ops.CorrFrame,2);
            badi = getOutliers(ops);
            ops1{i}.badframes(badi) = true;
            
            ind = repmat(~ops1{i}.badframes',1,2,nplanes) & ...
                repmat(permute(ops1{i}.usedPlanes, [1 3 2]),1,2,1);
            ds = ops1{i}.DS;
            ds = ds(ind); % [time x [x,y] x nplanes]
            ds = reshape(permute(ds, [1 3 2]), [], 2);
        else
            % determine bad frames
            badi                    = getOutliers(ops1{i});
            ops1{i}.badframes(badi) = true;
            
            ds = ops1{i}.DS(~ops1{i}.badframes,:);
        end
        
        minDs = min(ds, [], 1);
        maxDs = max(ds, [], 1);
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
end
    
savepath = sprintf('%s/', ops.ResultsSavePath);
if ~exist(savepath, 'dir')
    mkdir(savepath)
end
save(sprintf('%s/regops_%s_%s.mat', ops.ResultsSavePath, ...
    ops.mouse_name, ops.date),  'ops1')


%save(sprintf('%s/F_%s_%s_plane%d.mat', ops.ResultsSavePath, ...
 %   ops.mouse_name, ops.date, iplane), 'ops')

%%
