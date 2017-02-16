% computes registration offsets for data split into blocks
% loops over blocks and returns offsets dsall
function [dsall,ops1] = GetBlockOffsets(data, j, iplane0, ops, ops1)
            
nblocks = ops.numBlocks(1)*ops.numBlocks(2);
dsall = zeros(size(data,3), 2, nblocks);
for i = 1:ops.numPlanes
    ifr0 = iplane0(ops.planesToProcess(i));
    indframes = ifr0:ops.nplanes:size(data,3);
    ds = zeros(numel(indframes), 2, nblocks,'double');
    Corr = zeros(numel(indframes), nblocks,'double');
    for ib = 1:nblocks
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
    