function redoRegistration(ops, ichannel)

% Cannot be used with multiple FOVs, and ignores BiDiPhase

alignAcrossPlanes  = getOr(ops, {'alignAcrossPlanes'}, false);
interpolateAcrossPlanes = getOr(ops, {'interpolateAcrossPlanes'}, false);
planesToInterpolate = getOr(ops, {'planesToInterpolate'}, 1:nplanes);

if nargin < 2 || isempty(ichannel)
    if isfield(ops, 'gchannel') && ~isempty(ops.gchannel)
        ichannel = ops.gchannel;
    else
        ichannel = 1;
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

iplane = ops.iplane;
dsall = ops.DS;
T = size(dsAll,1);
Dims = size(dsAll,2);
usedPlanes = [];
if isfield(ops, 'alignedToPlanes')
    if ops.planeInterpolated == 0
        subscripts = cat(3, repmat((1:T)',1,Dims), ...
            repmat(1:Dims,T,1), repmat(ops.alignedToPlanes,1,Dims));
        subscripts = reshape(subscripts, [], 3);
        linearInd = sub2ind(size(dsall), subscripts(:,1), subscripts(:,2), ...
            subscripts(:,3));
        values = dsall(linearInd);
        values = reshape(values, T, Dims);
        dsall = values;
    else
        usedPlanes = find(sum(ops.usedPlanes,1) > 0);
        shifts = iplane - ops.basisPlanes;
        ds = NaN(T*usedPlanes, Dims);
        for i = 1:length(usedPlanes)
            subscripts = cat(3, repmat((1:T)',1,Dims), ...
                repmat(1:Dims,T,1), repmat(usedPlanes(i) + shifts,1,Dims));
            subscripts = reshape(subscripts, [], 3);
            subscripts(subscripts(:,3)<planesToInterpolate(1), 3) = ...
                planesToInterpolate(1);
            subscripts(subscripts(:,3)>planesToInterpolate(end), 3) = ...
                planesToInterpolate(end);
            linearInd = sub2ind(size(dsall), subscripts(:,1), subscripts(:,2), ...
                subscripts(:,3));
            values = dsall(linearInd);
            values = reshape(values, T, Dims);
            indframes = i:length(usedPlanes):(T*length(usedPlanes));
            ds(indframes,:) = values;
        end
        dsall = ds;
    end
end

t2 = 0;
% if two consecutive files have as many bytes, they have as many frames
nbytes = 0;
for k = 1:length(fs)
    ds = dsall(sum(ops.Nframes((1:length(fs))<k))*max(1,length(usedPlanes)) + ...
        (1:ops.Nframes(k)),:);
%     ds = ds(iplane0:ops.nplanes:end,:);
    dataPrev = zeros(ops.Ly, ops.Lx, nplanes, 'int16');
    iplane0 = 1:1:ops.nplanes;
    t1 = 0;
    for j = 1:length(fs{k})
        if abs(nbytes - fs{k}(j).bytes)>1e3
            nbytes = fs{k}(j).bytes;
            nFr = nFramesTiff(fs{k}(j).name);
        end
        iplane0 = mod(iplane0-1, nplanes) + 1;
        
        if isempty(usedPlanes)
            ichanset = [iplane0(iplane); nFr; nchannels*ops.nplanes];
        else
            ichanset = [ichannel; nFr; nchannels];
        end
        data = loadFramesBuff(fs{k}(j).name, ichanset(1), ichanset(2), ...
            ichanset(3));
        if ~isempty(usedPlanes)
            [~,a] = sort(iplane0(usedPlanes));
            planes = repmat(usedPlanes(a)', 1, ceil(size(data,3)/nplanes));
            planes = planes(:);
            ind = bsxfun(@plus, iplane0(usedPlanes), ...
                (0:nplanes:ceil(size(data,3)/nplanes)));
            ind = ind(:);
            planes(ind>size(data,3)) = [];
            ind(ind>size(data,3)) = [];
            data = data(:,:,ind);
        end
        ix0 = 0;
        Nbatch = 1000;
        dreg = zeros(size(data), class(data));
        while ix0<size(data,3)
            indxr = ix0 + (1:Nbatch);
            indxr(indxr>size(data,3)) = [];
            if iscell(splitBlocks)
                dreg(:, :, indxr) = blockRegisterMovie(data(:, :, indxr), ...
                    xyMask, ds(t1 + indxr,:,:));
            else
                dreg(:, :, indxr)        = ...
                    register_movie(data(:, :, indxr), ops, ds(t1 + indxr,:));
            end
            ix0 = ix0 + Nbatch;
        end
        if ~isempty(usedPlanes)
            ind1 = find(planes == iplane);
            bases = ops.basisPlanes(ceil((t2+t1)/length(usedPlanes) + ...
                (1:length(ind1))));
            uniqueBases = unique(bases);
            dwrite = zeros(size(dreg,1), size(dreg,2), length(ind1));
            dataNext = [];
            for b = uniqueBases
                planesInv = find(~isnan(ops.profiles(:,b)))';
                diffs = planesInv - iplane;
                for pl = 1:length(planesInv)
                    ind2 = ind1 + diffs(pl);
                    ind2(bases ~= b) = [];
                    weight = ops.profiles(planesInv(pl),b);
                    if ind2(1)<1 % frame is part of previously loaded tiff file (can happen if frames per tiff are not mupltiple of planes)
                        dwrite(:,:,1) = dwrite(:,:,1) + ...
                            weight .* double(dataPrev(:,:,end-ind2(1)));
                        ind2(1) = [];
                    end
                    if ind2(end)>size(dreg,3) % frame is part of upcoming tiff file
                        if j<length(fs{k})
                            if isempty(dataNext) % load first few frames of next tiff and register them
                                ichanset = [ichannel; nplanes; ops.nchannels];
                                data = loadFramesBuff(fs{k}(j+1).name, ichanset(1), ichanset(2), ichanset(3));
                                iplane1 = iplane0 - nFr/nchannels;
                                iplane1 = mod(iplane1-1, nplanes) + 1;
                                data = data(:,:,iplane1(usedPlanes));
                                dataNext = register_movie(data, ops, ...
                                    ds(t1+t2+length(planes)+ (1:length(usedPlanes)),:));
                            end
                            dwrite(:,:,end) = dwrite(:,:,end) + weight .*...
                                double(dataNext(:,:,ind2(end)-size(dreg,3)));
                        end
                        ind2(end) = [];
                    end
                    dwrite(:,:,ceil(ind2/length(usedPlanes))) = dwrite(:,:,ceil(ind2/length(usedPlanes))) + ...
                        weight .* double(dreg(:,:,ind2));
                end
            end
            dataPrev = dreg(:,:,end-length(usedPlanes)+1:end);
            fwrite(fid, dwrite, class(data));
            t1 = t1 + length(planes);
        else
            fwrite(fid, dreg, class(data));
            t1 = t1 + size(data,3);
        end
        iplane0 = iplane0 - nFr/nchannels;
    end
    t2 = t2 + ops.Nframes(k)*length(usedPlanes);
end

fclose(fid);