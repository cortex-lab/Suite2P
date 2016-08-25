function redoRegistration(ops, ichannel)

if nargin < 2 || isempty(ichannel)
    if isfield(ops, 'gchannel') && ~isempty(ops.gchannel)
        ichannel = ops.gchannel;
    else ichannel = 1;
    end
end

fs = ops.fsroot;
nchannels = ops.nchannels;

regdir = fileparts(ops.RegFile);
if ~exist(regdir, 'dir')
    mkdir(regdir);
end

splitBlocks    = getOr(ops, {'splitBlocks'}, 'none');
Ly = ops.Ly;
Lx = ops.Lx;
if iscell(splitBlocks)
    numBlocks = length(ops.yBL);
    xyMask = zeros(Ly, Lx, numBlocks, 'single');
    for i = 1:numBlocks
        msk = zeros(Ly, Lx, 'single');
        msk(ops.yBL{i},ops.xBL{i}) = 1;
        sT = ops.splitTaper;
        msk = my_conv(my_conv(msk, sT)',sT)';
        xyMask(:,:,i) = msk;
    end
    xyMask = xyMask./repmat(sum(xyMask, 3), 1, 1, numBlocks);
    xyMask = reshape(xyMask, Ly*Lx, numBlocks);
end

% open bin file for writing
fid = fopen(ops.RegFile, 'w');

iplane0 = ops.iplane;
dsall = ops.DS;
for k = 1:length(fs)
    ds = dsall(sum(ops.Nframes((1:length(fs))<k)) + (1:ops.Nframes(k)),:,:);
%     ds = ds(iplane0:ops.nplanes:end,:);
    frCount = 0;
    frCount2 = 0;
    for j = 1:length(fs{k})
        nFr = nFrames(fs{k}(j).name);
        frameGroups = mod(frCount+(1:nchannels*ops.nplanes), ...
            nchannels*ops.nplanes);
        frameGroups(frameGroups == 0) = nchannels*ops.nplanes;
        firstFrame = find(frameGroups == nchannels*(iplane0-1)+ichannel);
        ichanset = [firstFrame; nFr; nchannels*ops.nplanes];
        data = loadFramesBuff(fs{k}(j).name, ichanset(1), ichanset(2), ...
            ichanset(3));
        ix0 = 0;
        Nbatch = 1000;
        dreg = zeros(size(data), class(data));
        while ix0<size(data,3)
            indxr = ix0 + (1:Nbatch);
            indxr(indxr>size(data,3)) = [];
            if iscell(splitBlocks)
                dreg(:, :, indxr) = blockRegisterMovie(data(:, :, indxr), ...
                    xyMask, ds(frCount2 + indxr,:,:));
            else
                dreg(:, :, indxr)        = ...
                    register_movie(data(:, :, indxr), ops, ds(frCount2 + indxr,:));
            end
            ix0 = ix0 + Nbatch;
        end
        fwrite(fid, dreg, class(data));
        frCount = frCount + nFr;
        frCount2 = frCount2 + size(data,3);
    end
end

fclose(fid);