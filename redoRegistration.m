function redoRegistration(ops)

if isfield(ops, 'gchannel') && ~isempty(ops.gchannel); ichannel    = ops.gchannel;
else ichannel = 1; end

fs = ops.fsroot;
nchannels = ops.nchannels;

regdir = fileparts(ops.RegFile);
if ~exist(regdir, 'dir')
    mkdir(regdir);
end

% open bin file for writing
fid = fopen(ops.RegFile, 'w');

iplane0 = ops.iplane;
dsall = ops.DS;
for k = 1:length(fs)
    ds = dsall(sum(ops.Nframes((1:length(fs))<k)) + (1:ops.Nframes(k)),:);
%     ds = ds(iplane0:ops.nplanes:end,:);
    frCount = 0;
    frCount2 = 0;
    for j = 1:length(fs{k})
        nFr = img.nFrames(fs{k}(j).name);
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
            dreg(:, :, indxr)        = ...
                register_movie(data(:, :, indxr), ops, ds(frCount2 + indxr,:));
            ix0 = ix0 + Nbatch;
        end
        fwrite(fid, dreg, class(data));
        frCount = frCount + nFr;
        frCount2 = frCount2 + size(data,3);
    end
end

fclose(fid);