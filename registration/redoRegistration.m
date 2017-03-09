function redoRegistration(ops1, ichannel, iplanes)

% Cannot be used with multiple FOVs, and ignores BiDiPhase

if nargin < 2 || isempty(ichannel)
    if isfield(ops, 'gchannel') && ~isempty(ops.gchannel)
        ichannel = ops.gchannel;
    else
        ichannel = 1;
    end
end

fs = ops1{1}.fsroot;
nchannels = ops1{1}.nchannels;
nplanes = ops1{1}.nplanes;

regdir = fileparts(ops1{1}.RegFile);
if ~exist(regdir, 'dir')
    mkdir(regdir);
end

splitBlocks    = getOr(ops1{1}, {'splitBlocks'}, 'none');
Ly = ops1{1}.Ly;
Lx = ops1{1}.Lx;
if iscell(splitBlocks)
    numBlocks = length(ops1{1}.yBL);
    xyMask = zeros(Ly, Lx, numBlocks, 'single');
    for i = 1:numBlocks
        msk = zeros(Ly, Lx, 'single');
        msk(ops1{1}.yBL{i},ops1{1}.xBL{i}) = 1;
        sT = ops1{1}.splitTaper;
        msk = my_conv(my_conv(msk, sT)',sT)';
        xyMask(:,:,i) = msk;
    end
    xyMask = xyMask./repmat(sum(xyMask, 3), 1, 1, numBlocks);
    xyMask = reshape(xyMask, Ly*Lx, numBlocks);
end

% open bin file for writing
fid = cell(1, length(iplanes));
for i = 1:length(iplanes)
    fid{i} = fopen(ops1{iplanes(i)}.RegFile, 'w');
end

dsall = NaN(size(ops1{1}.DS,1), 2, nplanes);
for i = 1:nplanes
    t = 0;
    for k = 1:length(ops1{i}.Nframes)
        dsall(t + (1:ops1{i}.Nframes(k)),:,i) = ...
            ops1{i}.DS(sum(ops1{i}.Nframes(1:k-1)) + (1:ops1{i}.Nframes(k)),:);
        t = t + ops1{1}.Nframes(k);
    end
end
dsall = reshape(permute(dsall, [3 1 2]), [], 2);

t2 = 0;
% if two consecutive files have as many bytes, they have as many frames
nbytes = 0;
for k = 1:length(fs)
    dataPrev = zeros(Ly, Lx, nplanes, 'int16');
    iplane0 = 1:1:nplanes;
    t1 = 0;
    for j = 1:length(fs{k})
        if abs(nbytes - fs{k}(j).bytes)>1e3
            nbytes = fs{k}(j).bytes;
            nFr = nFramesTiff(fs{k}(j).name);
        end
        iplane0 = mod(iplane0-1, nplanes) + 1;
        
        ichanset = [ichannel; nFr; nchannels];
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
                    xyMask, dsall(t1+t2 + indxr,:,:));
            else
                dreg(:, :, indxr)        = ...
                    register_movie(data(:, :, indxr), ops1{1}, dsall(t1+t2 + indxr,:));
            end
            ix0 = ix0 + Nbatch;
        end
        
        dataNext = [];
        for i = 1:length(iplanes)
            if isfield(ops1{iplanes(i)}, 'planeInterpolated') && ...
                    ops1{iplanes(i)}.planeInterpolated == 1
                ind1 = iplane0(iplanes(i)) : nplanes : size(data,3);
                bases = ops1{iplanes(i)}.basisPlanes(ceil((t2+t1)/nplanes + ...
                    (1:length(ind1))));
                uniqueBases = unique(bases)';
                dwrite = zeros(size(dreg,1), size(dreg,2), length(ind1));
                for b = uniqueBases
                    planesInv = find(~isnan(ops1{iplanes(i)}.profiles(:,b)) & ...
                        ops1{iplanes(i)}.profiles(:,b)>0)';
                    diffs = planesInv - iplanes(i);
                    for pl = 1:length(planesInv)
                        ind2 = ind1 + diffs(pl);
                        ind2(bases ~= b) = [];
                        weight = ops1{iplanes(i)}.profiles(planesInv(pl),b);
                        if ind2(1)<1 % frame is part of previously loaded tiff file (can happen if frames per tiff are not mupltiple of planes)
                            dwrite(:,:,1) = dwrite(:,:,1) + ...
                                weight .* double(dataPrev(:,:,end-ind2(1)));
                            ind2(1) = [];
                        end
                        if ind2(end)>size(dreg,3) % frame is part of upcoming tiff file
                            if j<length(fs{k})
                                if isempty(dataNext) % load first few frames of next tiff and register them
                                    ichanset = [ichannel; nplanes; nchannels];
                                    data = loadFramesBuff(fs{k}(j+1).name, ...
                                        ichanset(1), ichanset(2), ichanset(3));
                                    dataNext = register_movie(data, ops1{iplanes(i)}, ...
                                        dsall(t1+t2+size(dreg,3)+ (1:nplanes),:));
                                end
                                dwrite(:,:,end) = dwrite(:,:,end) + weight .*...
                                    double(dataNext(:,:,ind2(end)-size(dreg,3)));
                            end
                            ind2(end) = [];
                        end
                        dwrite(:,:,ceil(ind2/nplanes)) = dwrite(:,:,ceil(ind2/nplanes)) + ...
                            weight .* double(dreg(:,:,ind2));
                    end
                end
                fwrite(fid{i}, dwrite, class(data));
            else
                fwrite(fid{i}, dreg, class(data));
            end
        end
        dataPrev = dreg(:,:,end-nplanes+1:end);
        t1 = t1 + size(dreg,3);
        iplane0 = iplane0 - nFr/nchannels;
    end
    t2 = t2 + ceil(t1/nplanes)*nplanes;
end

for i = 1:length(iplanes)
    fclose(fid{i});
end