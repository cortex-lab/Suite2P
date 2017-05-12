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

Ly = ops1{1}.Ly;
Lx = ops1{1}.Lx;
BiDiPhase = getOr(ops1{1}, {'BiDiPhase'}, 0);

% open bin file for writing
fid = cell(1, length(iplanes));
for i = 1:length(iplanes)
    fid{i} = fopen(ops1{iplanes(i)}.RegFile, 'w');
end

sz = size(ops1{1}.DS);
dsall = NaN([sz, nplanes]);
for i = 1:nplanes
    t = 0;
    for k = 1:length(ops1{i}.Nframes)
        sub = t + (1:ops1{i}.Nframes(k))';
        for d = 2:length(sz)
            sub = [repmat(sub, sz(d), 1), ...
                reshape(repmat(1:sz(d), length(sub), 1), [], 1)];
        end
        sub = [sub, ones(size(sub,1),1).*i];
        sub = mat2cell(sub, size(sub,1), ones(1, size(sub,2)));
        ind = sub2ind(size(dsall), sub{:});
        dsall(ind) = ops1{i}.DS(sum(ops1{i}.Nframes(1:k-1)) + ...
            (1:ops1{i}.Nframes(k)),:);
        t = t + ops1{1}.Nframes(k);
    end
end
dsall = reshape(permute(dsall, [length(sz)+1 1:length(sz)]), [nplanes*sz(1) sz(2:end)]);

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
        
        if abs(BiDiPhase) > 0
            data = ShiftBiDi(BiDiPhase, data, Ly, Lx);
        end
        
        if isfield(ops1{1}, 'numBlocks') && ~isempty(ops1{1}.numBlocks)
            dreg = BlockRegMovie(data, ops1{1}, ...
                dsall(t1+t2 + (1:size(data,3)),:,:), true(Ly,Lx));
        else
            dreg = register_movie(data, ops1{1}, dsall(t1+t2 + (1:size(data,3)),:));
        end
        
        dataNext = [];
        for i = 1:length(iplanes)
            ind1 = iplane0(iplanes(i)) : nplanes : size(data,3);
            if isfield(ops1{iplanes(i)}, 'planeInterpolated') && ...
                    ops1{iplanes(i)}.planeInterpolated == 1
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
                fwrite(fid{i}, dreg(:,:,ind1), class(data));
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