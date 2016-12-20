function [dsall,ops1] = GetBlockOffsets(data, j, iplane0, ops, ops1)
            
numBlocks = ops.numBlocks;
dsall = zeros(size(data,3), 2, ops.numBlocks);
for i = 1:ops.numPlanes
    ifr0 = iplane0(ops.planesToProcess(i));
    indframes = ifr0:ops.nplanes:size(data,3);
    ds = zeros(numel(indframes), 2, numBlocks,'double');
    Corr = zeros(numel(indframes), numBlocks,'double');
    for ib = 1:ops.numBlocks
        % collect ds
        ops1{i}.mimg = ops1{i}.mimgB{ib};
        [ds(:,:,ib), Corr(:,ib)]  = ...
            regoffKriging(data(ops1{i}.yBL{ib},ops1{i}.xBL{ib},...
            indframes), ops1{i}, 0);
    end
    if j==1
        ds(1,:,:) = 0;
    end
    dsall(indframes,:,:)  = ds;
    ops1{i}.DS          = cat(1, ops1{i}.DS, ds);
    ops1{i}.CorrFrame   = cat(1, ops1{i}.CorrFrame, Corr);
end
    