function [ops1, planesToInterp, xyValid] = registerBlocksAcrossPlanes( ...
    ops1, ops, fid, fidIntpol)

nplanes = getOr(ops, {'nplanes'}, 1);
planesToInterpolate = getOr(ops, {'planesToInterpolate'}, 1:nplanes);
ichannel = getOr(ops, {'gchannel'}, 1);
BiDiPhase = getOr(ops, {'BiDiPhase'}, 0);
RegFileBinLocation = getOr(ops, {'RegFileBinLocation'}, []);
interpolateAcrossPlanes = getOr(ops, {'interpolateAcrossPlanes'}, false);
fs = ops.fsroot;
numBlocks = prod(ops.numBlocks);

% in order to determine optimal shifts in z at each time point, collect 
% registration offsets and correlations of all frames with all target images
corrsAll = []; % [time x alignedPlanes (numPlanes) x 
               %  targetPlanes (length(planesToInterpolate)) x numBlocks]
dsAllPlanes = []; % [time x [x,y] x alignedPlanes (numPlanes) x 
                  %  targetPlanes (length(planesToInterpolate)) x numBlocks];
for i = planesToInterpolate
    % reshape ds to [frames x 2 x planes x numBlocks], and
    % Corr to [frames x planes x numBlocks]; (need to add NaNs at end
    % of each experiment if last frame was not recorded in last plane)
    ds = [];
    corrs = [];
    t = 0;
    for k = 1:length(ops1{i}.Nframes)
        d = ops1{i}.DS(t + (1:ops1{i}.Nframes(k)),:,:);
        c = ops1{i}.CorrFrame(t + (1:ops1{i}.Nframes(k)),:);
        framesToAdd = ceil(size(d,1)/nplanes)*nplanes - size(d,1);
        d = permute(reshape([d; NaN(framesToAdd,2,size(d,3))], ...
            nplanes, [], 2, numBlocks), [2 3 1 4]);
        c = permute(reshape([c; NaN(framesToAdd,size(c,2))], ...
            nplanes, [], numBlocks), [2 1 3]);
        ds = [ds; d];
        corrs = [corrs; c];
        t = t + ops1{i}.Nframes(k);
        fr = floor(ops1{i}.Nframes(k)/nplanes) + ...
            double(mod(ops1{i}.Nframes(k), nplanes) >= i);
        ops1{i}.Nframes(k) = fr;
    end
    ops1{i}.DS = ds; % [time x 2 x nplanes x numBlocks]
    ops1{i}.CorrFrame = corrs; % [time x nplanes x numBlocks]
    
    if ~isempty(dsAllPlanes) && size(ds,1)<size(dsAllPlanes,1)
        ds(end+1:size(dsAllPlanes,1),:,:,:) = NaN;
        corrs(end+1:size(corrsAll,1),:,:) = NaN;
    end
    dsAllPlanes = cat(4, dsAllPlanes, permute(ds, [1 2 3 5 4]));
    corrsAll = cat(3, corrsAll, permute(corrs, [1 2 4 3]));
end
for i = setdiff(2:nplanes, planesToInterpolate)
    ds = [];
    corrs = [];
    t = 0;
    for k = 1:length(ops1{i}.Nframes)
        d = NaN(ops1{1}.Nframes(k),2,numBlocks);
        d(1:ops1{i}.Nframes(k),:,:) = ops1{i}.DS(t + (1:ops1{i}.Nframes(k)),:,:);
        c = NaN(ops1{1}.Nframes(k),numBlocks);
        c(1:ops1{i}.Nframes(k),:) = ops1{i}.CorrFrame(t + (1:ops1{i}.Nframes(k)),:);
        ds = [ds; d];
        corrs = [corrs; c];
        t = t + ops1{i}.Nframes(k);
    end
    ops1{i}.DS = ds;
    ops1{i}.CorrFrame = corrs;
end

% find z-shifts in data; assume that whole stack is moving together
% (won't detect very fast z-shifts)
% only consider non-flyback planes here and average across all blocks of
% plane
corrsAll = mean(corrsAll(:,planesToInterpolate,:,:),4);
corrsNorm = bsxfun(@rdivide, corrsAll, max(corrsAll,[],2));
% for each possible shift in z, average alignment correlations across all
% planes (except those that are shifted outside of stack)
diags = -size(corrsAll,3)+2 : size(corrsAll,3)-2; % ignore maximum possible shift
elems = 2:size(corrsAll,3);
elems = [elems, flip(elems(1:end-1))];
diagonals = zeros(size(corrsAll,2),size(corrsAll,3),length(diags));
for d = 1:length(diags)
    diagonals(:,:,d) = diag(ones(1,elems(d)),diags(d));
end
diagonals = diagonals==1;
sh = zeros(size(corrsAll,1),length(diags));
for t = 1:size(corrsAll,1)
    t1 = squeeze(corrsNorm(t,:,:));
    for d = 1:length(diags)
        dgnl = diagonals(:,:,d);
        sh(t,d) = nansum(t1(dgnl))./sum(~isnan(t1(dgnl)));
    end
end
% for each frame, determine shift resulting in maximal alignment correlation
[corrs,inds] = max(sh,[],2);
shifts = diags(inds); % use shifts like this: bestTargetImage(t) = plane + shifts(t)
                      %                       bestPlane(t) = targetImage - shifts(t)
slowShifts = round(medfilt1(shifts,500));
ops2.CorrFrame = corrs;
outliers = getOutliers(ops2);
outliers = union(outliers, find(abs(shifts-slowShifts) > ...
    round(length(planesToInterpolate)/5)));
shifts(outliers) = slowShifts(outliers);

if getOr(ops, 'fig', 1) == 1
    figure
    plot(shifts)
    hold on
    yLimits = get(gca, 'YLim');
    plot(repmat(cumsum(ops1{1}.Nframes(1:end-1)),2,1), ...
        repmat(yLimits(:), 1, length(ops1{1}.Nframes)-1), 'r')
    xlabel('# Frames')
    ylabel('Shift in z direction')
end

% determine those planes that can be followed through the whole
% recording (while accounting for z-shifts)
indInterpolate = 1+max(shifts) : ...
    length(planesToInterpolate)+min(shifts);
if isempty(indInterpolate)
    disp('NOTE: no plane will be interpolated because range of shifts is too large')
end
planesToInterp = planesToInterpolate(indInterpolate); % target images
% of these planes will be used (at every time, there is one frame whose
% best matching target image is one of these target images)

% determine which planes are highly correlated to each other
% (similar neihgbouring planes); those will be averaged
corrsPadded = [NaN(size(corrsNorm)), corrsNorm, NaN(size(corrsNorm))];
for t = 1:size(corrsNorm,1)
    corrsPadded(t,:,:) = circshift(corrsPadded(t,:,:),shifts(t),2);
end
profiles = squeeze(nanmean(corrsPadded(:,size(corrsNorm,2)+ ...
    (1:size(corrsNorm,2)),:), 1));
ind = profiles>=0.8;
sngls = find(sum(ind,1) < 2);
for k = sngls
    [~,j] = sort(profiles(:,k),'descend');
    ind(j(2),k) = true;
end
profiles(~ind) = NaN;
profiles = bsxfun(@rdivide, profiles, nansum(profiles,1)); % [planesToInterpolate x planesToInterpolate]
if getOr(ops, 'fig', 1) == 1
    figure
    imagesc(profiles, [0 max(profiles(:))])
    xlabel('Target images')
    ylabel('Aligned planes')
    title('Z-profiles')
end

% for each frame, select registration offset: align to best matching
% target image
T = size(dsAllPlanes,1);
Dims = size(dsAllPlanes,2);
dsall = NaN(T*nplanes, 2, numBlocks);
% (1) collect dsall for flyback planes (not in planesToInterpolate)
for i = setdiff(1:nplanes, planesToInterpolate)
    indframes = ops.planesToProcess(i):nplanes:(T*nplanes);
    ds = ops1{i}.DS;
    indframes(size(ds,1)+1:end) = [];
    dsall(indframes,:,:) = ds;
    
    ops1{i}.planeInterpolated = 0;
    up = false(size(ds,1), nplanes);
    up(:,i) = true;
    ops1{i}.usedPlanes = up;
end
% (2) collect dsall for all planes that will be interpolated
for i = 1:length(planesToInterpolate)
    indframes = planesToInterpolate(i):nplanes:(T*nplanes);
    ds = squeeze(dsAllPlanes(:,:,planesToInterpolate(i),:,:));
    subscripts = cat(3, repmat((1:T)',1,Dims), ...
        repmat(1:Dims,T,1), repmat((i + shifts)',1,Dims));
    subscripts = reshape(subscripts, [], 3);
    subscripts(subscripts(:,3)<1, 3) = 1;
    subscripts(subscripts(:,3)>length(planesToInterpolate), 3) = ...
        length(planesToInterpolate);
    linearInd = sub2ind(size(ds(:,:,:,1)), subscripts(:,1), subscripts(:,2), ...
        subscripts(:,3));
    
    if ismember(i, indInterpolate)
        prof = eye(nplanes);
        prof(prof==0) = NaN;
        prof(planesToInterpolate, planesToInterpolate) = profiles;
        up = planesToInterpolate(i) - shifts';
        up = prof(:,up)';
        up = ~isnan(up);
    end
    for j = 1:numBlocks
        d = ds(:,:,:,j);
        values = d(linearInd);
        values = reshape(values, T, Dims);
        dsall(indframes,:,j) = values;
    end
    ops1{planesToInterpolate(i)}.DS_allPlanes = dsAllPlanes;
    ops1{planesToInterpolate(i)}.CorrFrame_allPlanes = corrsAll;
    values = dsall(indframes,:,:);
    nanInds = isnan(values(:,1,1));
    values(nanInds,:,:) = [];
    ops1{planesToInterpolate(i)}.DS = values;
    
    if ismember(i, indInterpolate)
        ops1{planesToInterpolate(i)}.planeInterpolated = 1;
        ops1{planesToInterpolate(i)}.usedPlanes = up;
        ops1{planesToInterpolate(i)}.basisPlanes = planesToInterpolate(i) - shifts';
        ops1{planesToInterpolate(i)}.basisPlanes(nanInds) = [];
        ops1{planesToInterpolate(i)}.profiles = prof;
    else
        ops1{planesToInterpolate(i)}.planeInterpolated = 0;
        up = false(size(ds,1), nplanes);
        up(:,planesToInterpolate(i)) = true;
        ops1{planesToInterpolate(i)}.usedPlanes = up;
    end
    ap = planesToInterpolate(i) + shifts';
    ap(ap<planesToInterpolate(1)) = planesToInterpolate(1);
    ap(ap>planesToInterpolate(end)) = planesToInterpolate(end);
    ops1{planesToInterpolate(i)}.alignedToPlanes = ap;
    
    ops =  ops1{planesToInterpolate(i)};
    ops.CorrFrame(repmat(~ops.usedPlanes,1,1,numBlocks)) = NaN;
    ops.CorrFrame = nanmean(nanmean(ops.CorrFrame,2),3); % average correlations across all used planes and all blocks
    ops.CorrFrame(nanInds) = [];
    ops1{planesToInterpolate(i)}.CorrFrame = ops.CorrFrame;
    ops1{planesToInterpolate(i)}.usedPlanes(nanInds,:) = [];
    ops1{planesToInterpolate(i)}.alignedToPlanes(nanInds) = [];
end

% align frames and write movie
t2 = 0;
% if two consecutive files have as many bytes, they have as many frames
nbytes = 0;
xyValid = true(ops.Ly, ops.Lx);
for k = 1:length(fs)
    dataPrev = zeros(ops.Ly, ops.Lx, nplanes, 'int16');
    iplane0 = 1:1:ops.nplanes;
    t1 = 0;
    for j = 1:length(fs{k})
        if abs(nbytes - fs{k}(j).bytes)>1e3
            nbytes = fs{k}(j).bytes;
            nFr = nFramesTiff(fs{k}(j).name);
        end
        iplane0 = mod(iplane0-1, nplanes) + 1;
        
        if mod(nFr, ops.nchannels) ~= 0
            fprintf('  WARNING: number of frames in tiff (%d) is NOT a multiple of number of channels!\n', j);
        end
        % load frames of all planes
        ichanset = [ichannel; nFr; ops.nchannels];
        data = loadFramesBuff(fs{k}(j).name, ichanset(1), ichanset(2), ichanset(3));
        if BiDiPhase
            yrange = 2:2:Ly;
            if BiDiPhase>0
                data(yrange, (1+BiDiPhase):Lx,:) = data(yrange, 1:(Lx-BiDiPhase),:);
            else
                data(yrange, 1:Lx+BiDiPhase,:) = data(yrange, 1-BiDiPhase:Lx,:);
            end
        end
        
        % align frames according to determined registration offsets
        [dreg, xyValid] = BlockRegMovie(data, ops, ...
            dsall(t1+t2+(1:size(data,3)),:,:), xyValid);
        
        % for each plane, write aligned movie to bin file+
        for i = 1:nplanes
            ifr0 = iplane0(ops.planesToProcess(i));
            indframes = ifr0:nplanes:size(data,3);
            dwrite = dreg(:,:,indframes);
            fwrite(fid{i}, dwrite, class(data));
            
            ops1{i}.mimg1 = ops1{i}.mimg1 + sum(dwrite,3);
        end
        
        if interpolateAcrossPlanes && ~isempty(RegFileBinLocation)
            % use frames of plane that match best to target image of
            % current plane and average across similar planes
            dataNext = [];
            for i = indInterpolate
                ind1 = iplane0(planesToInterpolate(i)) : nplanes : size(data,3);
                sh = shifts(ceil((t1+t2)/nplanes + (1:length(ind1))));
                uniqueShifts = unique(sh);
                dwrite = zeros(size(dreg,1),size(dreg,2),length(ind1));
                for s = uniqueShifts
                    planesInv = find(~isnan(profiles(:,i-s)))';
                    diffs = planesInv - i;
                    for pl = 1:length(planesInv)
                        ind2 = ind1 + diffs(pl); % smallest possible value: -nplanes+1,
                        % largest possible values: size(dreg,3)+nplanes
                        ind2(sh~=s) = [];
                        weight = profiles(planesInv(pl),i-s);
                        if ind2(1)<1 % frame is part of previously loaded tiff file (can happen if frames per tiff are not mupltiple of planes)
                            dwrite(:,:,1) = dwrite(:,:,1) + ...
                                weight .* double(dataPrev(:,:,end-ind2(1)));
                            ind2(1) = [];
                        end
                        if ind2(end)>size(dreg,3) % frame is part of upcoming tiff file
                            if j<length(fs{k})
                                if isempty(dataNext) % load first few frames of next tiff and register them
                                    ichanset = [ichannel; nplanes; ops.nchannels];
                                    data = loadFramesBuff(fs{k}(j+1).name, ...
                                        ichanset(1), ichanset(2), ichanset(3));
                                    if BiDiPhase
                                        yrange = 2:2:Ly;
                                        if BiDiPhase>0
                                            data(yrange, (1+BiDiPhase):Lx,:) = data(yrange, 1:(Lx-BiDiPhase),:);
                                        else
                                            data(yrange, 1:Lx+BiDiPhase,:) = data(yrange, 1-BiDiPhase:Lx,:);
                                        end
                                    end
                                    dataNext = register_movie(data, ops, ...
                                        dsall((t1+t2) + size(dreg,3) ...
                                        + (1:nplanes),:));
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
                fwrite(fidIntpol{ops.planesToProcess(planesToInterpolate(i))}, dwrite, class(data));
            end
        end
        dataPrev = dreg(:,:,end-nplanes+1:end);
        t1 = t1 + size(dreg,3);
        iplane0 = iplane0 - nFr/ops.nchannels;
    end
    t2 = t2 + ceil(t1/nplanes)*nplanes;
end