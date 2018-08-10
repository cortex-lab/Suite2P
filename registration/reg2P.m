function ops1 = reg2P(sm)

ops = sm.options;

%%
ops.doRegistration  = getOr(ops, {'doRegistration'}, 1); % register tiffs?

if ops.doRegistration
    disp('running registration');
else
    disp('skipping registration, but assembling binary file');
end

ops = buildRegOps(ops);
planesToProcess = ops.planesToProcess;
% --- number of planes in recording and number of channels --- %
nplanes        = ops.nplanes;
numPlanes          = length(planesToProcess); % # of planes to process
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

% Save options back to storage manager
sm.options = ops;

BiDiPhase              = ops.BiDiPhase;

% fs = ops.fsroot;

%% find the mean frame after aligning a random subset

% check if there are tiffs in directory
fileName = sm.getFirstFile();
try
    IMG = loadFramesBuff(fileName, 1, 1, 1);
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
    IMG = GetRandFrames(sm, ops);
    if isempty(IMG)
        error('ERROR: There are too few frames per plane for processing!');
    end

    % from random frames get scaling factor (if uint16, scale by 2)
    if ~isfield(ops, 'scaleTiff')
        if max(IMG(:)) > 2^15
            ops.scaleTiff = 2;
        else
            ops.scaleTiff = 1;
        end
    end

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
    ops1 = cell(numPlanes, size(xFOVs, 2));
    for i = 1:numPlanes
        for l = 1:size(xFOVs, 2)
            if ops.nonrigid && ~(ops.alignAcrossPlanes || ops.interpolateAcrossPlanes)
                ops1{i} = nonrigidAlignIterative(squeeze(IMG(:, :, planesToProcess(i), :)), ops);
            else
                ops1{i, l} = alignIterative(single(squeeze(IMG(yFOVs(:, l), xFOVs(:, l), ...
                    planesToProcess(i), :))), ops);
            end
        end
        fprintf('target image acquired for plane %d\n', planesToProcess(i));
    end

    if ops.alignTargetImages % align target images of all planes to each other
        % (reasonable only if interplane distance during imaging was small)
        newTargets = getImgOverPlanes(ops1);
        for i = 1:length(ops.planesToInterpolate)
            for l = 1:size(xFOVs, 2)
                ops1{ops.planesToInterpolate(i), l}.mimg = newTargets(:, :, i, l);
            end
        end
    end

    % display target image
    if ops.fig && ops.showTargetRegistration
        PlotRegMean(ops1, ops);
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
    for j = 1:size(xFOVs, 2)
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

sm.options = ops;
sm.resetCurrentEntry();

%% open files for registration
[ops1, fid, fidRED, fidIntpol] = openBinFiles(sm, ops1, red_binary);

%%
tic
% compute registration offsets and align using offsets
% if two consecutive files have as many bytes, they have as many frames
nbytes = 0;

xyValid = true(Ly, Lx);
sm.resetCurrentEntry();
lastPartition = -1;
fileCounter = 0;
tempTiffPath = sm.getOr('temp_tiff');

while sm.hasData()
    [fileName, fileInfo] = sm.getFile();
    fileCounter = fileCounter + 1;

    if fileInfo.partition ~= lastPartition
        iplane0 = 1:1:ops.nplanes; % identity of planes for first frames in tiff file
                                   % reset for each new partition

        nchannels_expt = sm.getNumChannels();
        if nchannels_expt > 1
            red_mean_expt = red_mean;
        else
            red_mean_expt = 0;
        end

        if red_align
            reg_channel = rchannel;
        else
            reg_channel = ichannel;
        end

        % initialize frame count
        for i = 1:numel(ops1)
            ops1{i}.Nframes(fileInfo.partition)     = 0;
        end
        lastPartition = fileInfo.partition;
    end

    % only compute number of frames if size of tiff is different
    % from previous tiff
    if abs(nbytes - fileInfo.bytes) > 1e3
        nbytes = fileInfo.bytes;
        nFr = nFramesTiff(fileName);
    end
    if mod(nFr, nchannels_expt) ~= 0
        fprintf('  WARNING: number of frames in tiff (%d) is NOT a multiple of number of channels!\n', nchannels_expt);
    end
    iplane0 = mod(iplane0-1, numPlanes) + 1;
    % only load frames of registration channel
    data = loadFramesBuff(fileName, reg_channel, nFr, nchannels_expt, tempTiffPath);
    if ~isempty(tempTiffPath)
        tiff_file = tempTiffPath;
    else
        tiff_file = fileName;
    end
    if abs(BiDiPhase) > 0
        data = ShiftBiDi(BiDiPhase, data, Ly, Lx);
    end

    % get the registration offsets for each frame
    if ops.doRegistration
        if ops.nonrigid
            [dsall, ops1] = nonrigidOffsets(data, lastPartition, fileCounter, iplane0, ops, ops1);
        else
            [dsall, ops1] = rigidOffsets(data, lastPartition, fileCounter, iplane0, ops, ops1);
        end
    end

    if ops.alignAcrossPlanes
        if rem(fileCounter, 5) == 1
            fprintf('Set %d, tiff %d done in time %2.2f \n', lastPartition, fileCounter, toc);
        end
        iplane0 = iplane0 - nFr/nchannels_expt;
        continue;
    end

    if ops.doRegistration
        % if aligning by the red channel, data needs to be reloaded as the green channel
        if red_align
            data_red = data;
            data = loadFramesBuff(tiff_file, ichannel, nFr, nchannels_expt, tempTiffPath);
            if abs(BiDiPhase) > 0
                data = ShiftBiDi(BiDiPhase, data, Ly, Lx);
            end
        end
        % load red channel for red binary file if not already loaded for red alignment
        if red_mean_expt && ~red_align
            data_red = loadFramesBuff(tiff_file, rchannel, nFr, nchannels_expt, tempTiffPath);
            if abs(BiDiPhase) > 0
                data_red = ShiftBiDi(BiDiPhase, data_red, Ly, Lx);
            end
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

    k = fileInfo.partition;

    % write dreg to bin file+
    for i = 1:numPlanes
        ifr0 = iplane0(planesToProcess(i));
        indframes = ifr0:nplanes:size(data, 3);
        for l = 1:size(xFOVs, 2)
            dwrite = dreg(yFOVs(:, l), xFOVs(:, l), indframes);
            dwrite = int16(dwrite / ops.scaleTiff);
            fwrite(fid{i, l}, dwrite, 'int16');
            ops1{i, l}.mimg1 = ops1{i, l}.mimg1 + sum(dwrite, 3);

            if red_mean_expt || red_align
                dwrite = dreg2(yFOVs(:, l), xFOVs(:, l), indframes);
                if ops.scaleTiff > 1
                    dwrite = int16(dwrite / ops.scaleTiff);
                end
                ops1{i, l}.mimgRED = ops1{i, l}.mimgRED + sum(dwrite, 3);
            end
            if red_binary
                fwrite(fidRED{i, l}, dwrite, 'int16');
            end
        end
    end

    if rem(fileCounter, 5) == 1
        fprintf('Set %d, tiff %d done in time %2.2f \n', lastPartition, fileCounter, toc);
    end

    iplane0 = iplane0 - nFr/nchannels_expt;
end

sm.options = ops;

if ops.alignAcrossPlanes && ops.doRegistration % align each frame with the best matching target image
    if ops.nonrigid
        [ops1, planesToInterp] = registerBlocksAcrossPlanes(ops1, sm, fid, fidIntpol);
    else
        [ops1, planesToInterp] = registerAcrossPlanes(ops1, ops, fid, fidIntpol);
    end
end

for i = 1:numel(ops1)
    ops1{i}.mimg1 = ops1{i}.mimg1/sum(ops1{i}.Nframes);
    if red_mean || red_align
        if ~isempty(ops.expts)
            red_expts = ismember(ops.expts, getOr(ops, 'expred', []));
        else
            red_expts = 1;
        end
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

% get mean of first and last frames in block (to check for drift)
if ~isempty(ops.nimgbegend) && ops.nimgbegend > 0
    for i = 1:numel(ops1)
        fid{i}           = fopen(ops1{i}.RegFile, 'r');
        ops1{i} = getBlockBegEnd(fid{i}, ops1{i}); % get mean of first and last frames in block (to check for drift)
        fclose(fid{i});
    end
end

% write registered tiffs to disk if ~isempty(sm.getOr('RegFileTiffLocation'))
if ~isempty(sm.getOr('RegFileTiffLocation'))
    for i = 1:numel(ops1)
        fid{i}           = fopen(ops1{i}.RegFile, 'r');
        write_reg_to_tiff(fid{i}, ops1{i}, i, 0, sm);
        fclose(fid{i});

        if red_binary
            fidRED{i}           = fopen(ops1{i}.RegFile2, 'r');
            write_reg_to_tiff(fidRED{i}, ops1{i}, i, red_binary, sm);
            fclose(fidRED{i});
        end
    end
end

regFileBinLocation = sm.getOr('RegFileBinLocation');
if ~isempty(regFileBinLocation)
    copyBinFile(ops1, sm);
end

if ops.interpolateAcrossPlanes && ~isempty(regFileBinLocation)
    filename = sm.registrationFileOnSlowStorage(i, 'interp');

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

sm.options = ops;
%options = ops1;
sm.saveRegistrationOptions(ops1);

%%
