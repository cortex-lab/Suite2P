function ops1 = reg2P(ops)
%%
ops.doRegistration  = getOr(ops, {'doRegistration'}, 1); % register tiffs?

if ops.doRegistration
    disp('running registration');
else
    disp('skipping registration, but assembling binary file');
end

ops = buildRegOps(ops);
% --- number of planes in recording and number of channels --- %
nplanes        = ops.nplanes;
numPlanes          = length(ops.planesToProcess); % # of planes to process
% which channel is the functional channel
ichannel           = getOr(ops, {'gchannel'}, 1);
% which channel is the non-functional channel
rchannel           = getOr(ops, {'rchannel'}, 2);

% --- red channel options ---%
red_align          = getOr(ops, {'AlignToRedChannel'}, 0); % register planes to red channel
red_binary          = getOr(ops, {'REDbinary'}, 0); % write red channel to a binary file
% extract mean red channel from blocks of recording with two channels
red_mean            = getOr(ops, {'redMeanImg'}, 0);
if red_binary
    red_mean = 1;
end
if red_mean
    disp('computing mean RED image if ~isempty(db.expred)');
end

BiDiPhase              = ops.BiDiPhase;

fs = ops.fsroot;

%% find the mean frame after aligning a random subset

% check if there are tiffs in directory
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
    % get frames for initial registration
    IMG = GetRandFrames(fs, ops);
    % compute phase shifts from bidirectional scanning
    if ops.dobidi
        ops.BiDiPhase = BiDiPhaseOffsets(IMG);
    end
    BiDiPhase = ops.BiDiPhase;
    fprintf('bi-directional scanning offset = %d pixels\n', BiDiPhase);
    if abs(BiDiPhase) > 0
        IMG = ShiftBiDi(BiDiPhase, IMG, Ly, Lx);
    end
    % makes blocks (number = numBlocks) and masks for smoothing registration offsets across blocks
    if ops.nonrigid
        ops = MakeBlocks(ops);
    end
    
    % for each plane: align chosen frames to average to generate target image
    ops1 = cell(numPlanes, size(xFOVs,2));
    for i = 1:numPlanes
        for l = 1:size(xFOVs,2)
            if ops.nonrigid && ~(ops.alignAcrossPlanes || ops.interpolateAcrossPlanes)
                ops1{i} = nonrigidAlignIterative(squeeze(IMG(:,:,ops.planesToProcess(i),:)), ops);
            else
                ops1{i,l} = alignIterative(single(squeeze(IMG(yFOVs(:,l),xFOVs(:,l),...
                    ops.planesToProcess(i),:))), ops);
            end
        end
        fprintf('target image acquired for plane %d\n', ops.planesToProcess(i));
    end
    
    if ops.alignTargetImages % align target images of all planes to each other
        % (reasonable only if interplane distance during imaging was small)
        newTargets = getImgOverPlanes(ops1, ops.planesToInterpolate);
        for i = 1:length(ops.planesToInterpolate)
            for l = 1:size(xFOVs,2)
                ops1{ops.planesToInterpolate(i),l}.mimg = newTargets(:,:,i,l);
            end
        end
    end
    
    % display target image
    if ops.fig && ops.showTargetRegistration
        PlotRegMean(ops1,ops);
        drawnow
    end
    clear IMG
    
else   % don't recompute mean image
    ops1 = cell(numPlanes, 1);
    for i = 1:numPlanes
        ops1{i} = ops;
        ops1{i}.mimg = zeros(Ly, Lx);
    end
end

%% initialize mean imgs and shifts
for i = 1:numPlanes
    for j = 1:size(xFOVs,2)
        if red_mean || red_align
            ops1{i,j}.mimgRED       = zeros(ops1{i,j}.Ly, ops1{i,j}.Lx);
        end
        ops1{i,j}.DS          = [];
        ops1{i,j}.CorrFrame   = [];
        ops1{i,j}.mimg1       = zeros(ops1{i,j}.Ly, ops1{i,j}.Lx);
        if ops.nonrigid
            ops1{i}.mimgB = cell(prod(ops.numBlocks),1);
            for ib = 1:ops.numBlocks(1)*ops.numBlocks(2)
                ops1{i}.mimgB{ib} = ops1{i}.mimg(ops1{i}.yBL{ib}, ops1{i}.xBL{ib});
            end
        end
    end
end

%% open files for registration

[ops1, fid, fidRED, fidIntpol] = openBinFiles(ops, ops1, red_binary);


%%
tic
% compute registration offsets and align using offsets
% if two consecutive files have as many bytes, they have as many frames
nbytes = 0;

xyValid = true(Ly,Lx);
for k = 1:length(fs)
    % in case different experiments have different numbers of channels
    if ismember(ops.expts(k), getOr(ops, 'expred', []))
        nchannels_expt = ops.nchannels_red;
        red_mean_expt  = red_mean;
    else
        nchannels_expt = ops.nchannels;
        red_mean_expt  = 0;
    end
    if red_align
        reg_channel = rchannel;
    else
        reg_channel = ichannel;
    end
        
    % initialize frame count
    for i = 1:numel(ops1)
        ops1{i}.Nframes(k)     = 0;
    end
    
    iplane0 = 1:1:ops.nplanes; % identity of planes for first frames in tiff file
    for j = 1:length(fs{k})
        % only compute number of frames if size of tiff is different
        % from previous tiff
        if abs(nbytes - fs{k}(j).bytes)>1e3
            nbytes = fs{k}(j).bytes;
            nFr = nFramesTiff(fs{k}(j).name);
        end
        if mod(nFr, nchannels_expt) ~= 0
            fprintf('  WARNING: number of frames in tiff (%d) is NOT a multiple of number of channels!\n', j);
        end
                
        iplane0 = mod(iplane0-1, numPlanes) + 1;
        % only load frames of registration channel
        data = loadFramesBuff(fs{k}(j).name, reg_channel, nFr, nchannels_expt, ops.temp_tiff);
        if ~isempty(ops.temp_tiff)
            tiff_file = ops.temp_tiff;
        else
            tiff_file = fs{k}(j).name;
        end
        if abs(BiDiPhase) > 0; data = ShiftBiDi(BiDiPhase, data, Ly, Lx); end
        
        % get the registration offsets for each frame
        if ops.doRegistration
            if ops.nonrigid
                [dsall, ops1] = nonrigidOffsets(data, j, iplane0, ops, ops1);
            else
                [dsall, ops1] = rigidOffsets(data, j, iplane0, ops, ops1);
            end
        end
        
        if ops.alignAcrossPlanes
            if rem(j,5)==1
                fprintf('Set %d, tiff %d done in time %2.2f \n', k, j, toc)
            end
            iplane0 = iplane0 - nFr/nchannels_expt;
            continue;
        end
        
        if ops.doRegistration
            % if aligning by the red channel, data needs to be reloaded as the green channel
            if red_align
                data_red = data;
                data = loadFramesBuff(tiff_file, ichannel, nFr, nchannels_expt, ops.temp_tiff);
                if abs(BiDiPhase) > 0; data = ShiftBiDi(BiDiPhase, data, Ly, Lx);  end
            end
            % load red channel for red binary file if not already loaded for red alignment
            if red_mean_expt && ~red_align
                data_red = loadFramesBuff(tiff_file, rchannel, nFr, nchannels_expt, ops.temp_tiff);
                if abs(BiDiPhase) > 0; data_red = ShiftBiDi(BiDiPhase, data_red, Ly, Lx); end
            end
            % align the frames according to the registration offsets
            if ops.nonrigid
                [dreg, xyValid] = nonrigidMovie(data, ops, dsall, xyValid);
            else
                dreg = rigidMovie(data, ops1, dsall, yFOVs, xFOVs);
            end
            
            if red_mean_expt || red_align
                if ops.nonrigid
                    [dreg2, xyValid] = nonrigidMovie(data_red, ops, dsall, xyValid);
                else
                    dreg2 = rigidMovie(data_red, ops1, dsall, yFOVs, xFOVs);
                end
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
                
                if red_mean_expt || red_align
                    dwrite = dreg2(yFOVs(:,l),xFOVs(:,l),indframes);
                    ops1{i,l}.mimgRED = ops1{i,l}.mimgRED + sum(dwrite,3);
                end
                if red_binary
                    fwrite(fidRED{i,l}, dwrite, class(data));
                end
            end
        end
        
        
        if rem(j,5)==1
            fprintf('Set %d, tiff %d done in time %2.2f \n', k, j, toc)
        end
        
        %keyboard;
        iplane0 = iplane0 - nFr/nchannels_expt;
    end
end

if ops.alignAcrossPlanes && ops.doRegistration % align each frame with the best matching target image
    if ops.nonrigid
        [ops1, planesToInterp] = registerBlocksAcrossPlanes(ops1, ops, fid, fidIntpol);
    else
        [ops1, planesToInterp] = registerAcrossPlanes(ops1, ops, fid, fidIntpol);
    end
end

for i = 1:numel(ops1)
    ops1{i}.mimg1 = ops1{i}.mimg1/sum(ops1{i}.Nframes);
    if red_mean || red_align
        red_expts = ismember(ops.expts, getOr(ops, 'expred', []));
        ops1{i}.mimgRED = ops1{i}.mimgRED/sum(ops1{i}.Nframes(red_expts));
        if red_binary
            fclose(fidRED{i});
        end
    end
    ops1{i}.badframes = false(1, size(ops1{i}.DS,1));
    if isfield(ops, 'badframes0') && ~isempty(ops.badframes0)
        ops1{i}.badframes(ops.badframes0) = true;
    end
    fclose(fid{i});
    if ops.interpolateAcrossPlanes
        fclose(fidIntpol{i});
    end
end

%%

% get mean of first and last frames in block (to check for drift)
if ~isempty(ops.nimgbegend) && ops.nimgbegend>0
    for i = 1:numel(ops1)
        fid{i}           = fopen(ops1{i}.RegFile, 'r');
        ops1{i} = getBlockBegEnd(fid{i}, ops1{i}); % get mean of first and last frames in block (to check for drift)
        fclose(fid{i});
    end
end


% write registered tiffs to disk if ~isempty(ops.RegFileTiffLocation)
if ~isempty(ops.RegFileTiffLocation)
    for i = 1:numel(ops1)
        fid{i}           = fopen(ops1{i}.RegFile, 'r');
        ops1{i} = write_reg_to_tiff(fid{i}, ops1{i}, i, 0);
        fclose(fid{i});
        
        if red_binary
            fidRED{i}           = fopen(ops1{i}.RegFile2, 'r');
            ops1{i} = write_reg_to_tiff(fidRED{i}, ops1{i}, i, 1);
            fclose(fidRED{i});
        end   
    end
end

% copy binary file if ~isempty(ops.RegFileBinLocation)
if ~isempty(ops.RegFileBinLocation)
    copyBinFile(ops1);
end

if ops.interpolateAcrossPlanes == 1 && ~isempty(RegFileBinLocation)
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

%%
% compute outlier frames and xrange, yrange of registered frames
if ops.doRegistration
    ops1 = getRangeNoOutliers(ops, ops1);
else
    for i = 1:numel(ops1)
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
