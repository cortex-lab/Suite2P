% computes registration offsets for data split into blocks
% loops over blocks and returns offsets dsall
function [dsall,ops1] = nonrigidOffsets(data, j, iplane0, ops, ops1)

nplanes = getOr(ops, {'nplanes'}, 1);
alignAcrossPlanes  = getOr(ops, {'alignAcrossPlanes'}, false);
planesToInterpolate = getOr(ops, {'planesToInterpolate'}, 1:nplanes);

nblocks = ops.numBlocks(1)*ops.numBlocks(2);
dsall = zeros(size(data,3), 2, nblocks);

for i = 1:numel(ops.planesToProcess)
    if alignAcrossPlanes && ismember(i, planesToInterpolate) % load frames of all planes
        indframes = 1:size(data,3);
    else
        ifr0 = iplane0(ops.planesToProcess(i));
        indframes = ifr0:ops.nplanes:size(data,3);
    end
    NT = numel(indframes);
    ds = zeros(NT, 2, nblocks,'double');
    Corr = zeros(NT, nblocks,'double');
    mimg0 = ops1{i}.mimg;
    for ib = 1:nblocks
        % collect ds
        ops1{i}.mimg = ops1{i}.mimgReg(ops1{i}.yBL{ib},ops1{i}.xBL{ib});
        if ops.kriging
            ups = 1;
            ops1{i}.SubPixel = 1;
            [ds(:,:,ib), Corr(:,ib), ccmat(:,:,ib)]  = ...
                regoffKriging(data(ops1{i}.yBL{ib},ops1{i}.xBL{ib},indframes),ops1{i}, 0);
        else
            [ds(:,:,ib), Corr(:,ib)]  = ...
                regoffLinear(data(ops1{i}.yBL{ib},ops1{i}.xBL{ib},indframes),ops1{i},0);
        end
    end
    ops1{i}.mimg = mimg0;
    
    %% smoothing for xblocks = 1
    
    ops1{i}.SubPixel = ops.SubPixel;
    subpixel = min(10,ops.SubPixel);
    %lcorr = (sqrt(size(ccmat,1))-1)/2;
    diffb = mean(ops1{i}.yBL{ops.numBlocks(2)+1}) - mean(ops1{i}.yBL{1});
    %keyboard;
    %%
    cc2 = reshape(ccmat,size(ccmat,1),NT,ops.numBlocks(2),ops.numBlocks(1));
    cc2 = permute(cc2, [1 3 2 4]);
    
    cc2 = cc2(:,:,:);
    sT = [1 1];
    cc2 = my_conv2(cc2, sT, [2 3]);
    cc2 = permute(reshape(cc2,size(cc2,1),ops.numBlocks(2),NT,ops.numBlocks(1)), [1 3 2 4]);
    cc2 = cc2(:,:,:);
    %%
    %keyboard;
    %%
    maxregshift = 5;
    %keyboard;
    % compute kernels for regression
    lpad     = 3;
    [x1,x2] = ndgrid([-lpad:lpad]);
    xt = [x1(:) x2(:)]';
    if ops.useGPU
        xt = gpuArray(single(xt));
    end

    sigL     = .85; % kernel width in pixels
    Kx = kernelD(xt,xt,sigL*[1;1]);
    linds = [-lpad:1/subpixel:lpad];
    [x1,x2] = ndgrid(linds);
    xg = [x1(:) x2(:)]';
    if ops.useGPU
        xg = gpuArray(single(xg));
    end
    Kg = kernelD(xg,xt,sigL*[1;1]);
    Kmat = Kg/Kx;
    dssmooth = [];
    for ib = 1:nblocks
        ly = numel(ops1{i}.yBL{ib});
        lx = numel(ops1{i}.xBL{ib});
        ccib = cc2(:,:,ib);
        ccib = reshape(ccib,ly,lx,NT);
        ccn  = -100*ones(size(ccib),'single');
        ccn(floor(ly/2)+1+[-maxregshift:maxregshift], floor(lx/2)+1+[-maxregshift:maxregshift], :) = ...
            ccib(floor(ly/2)+1+[-maxregshift:maxregshift], floor(lx/2)+1+[-maxregshift:maxregshift], :);
        [cmax,ii]   = max(reshape(ccn,[],NT),[],1);
        
        [iy, ix] = ind2sub([ly lx], ii);
        
        dl       = single(-lpad:1:lpad);
        if ops.useGPU
            dl   = gpuArray(dl);
            ccm = gpuArray.zeros(numel(dl), numel(dl), NT, 'single');
        else
            ccm = zeros(numel(dl), numel(dl), NT, 'single');
        end
        mxpt        = [iy(:) ix(:)];
        cccrop = [];
        
        for j = 1:NT
            % matrix +/- lpad surrounding max point
            cccrop(:, :, j) = ccib(mxpt(j,1)+dl, mxpt(j,2)+dl, j);
        end
        %
        cccrop = reshape(cccrop,[], NT);
        ccb         = Kmat * cccrop;
        
        % find max of grid
        [cx,ix]     = max(ccb, [], 1);
        [ix11,ix21] = ind2sub(numel(linds)*[1 1],ix);
        mdpt        = floor(numel(linds)/2)+1;
        dv0         = bsxfun(@minus, ([ix11' ix21'] - mdpt)/subpixel + mxpt, ...
            [floor(ly/2)+1 floor(lx/2)+1]);
        dssmooth(:,:,ib) = gather_try(dv0);
    end
  
    %%        
        
    if j==1
        dssmooth(1,:,:) = 0;
    end
    dsall(indframes,:,:)  = dssmooth;
    ops1{i}.DS          = cat(1, ops1{i}.DS, dssmooth);
    ops1{i}.CorrFrame   = cat(1, ops1{i}.CorrFrame, Corr);
end
