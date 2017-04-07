% registers frames using offsets dsall (rigid registration)
% loops over batches
function [dreg, xyValid] = BlockRegMovie(data, ops, dsall, xyValid)

ix0 = 0;
Nbatch = 1000;
dreg = zeros(size(data), 'like', data);
while ix0<size(data,3)
    indxr = ix0 + (1:Nbatch);
    indxr(indxr>size(data,3)) = [];
    
    [dreg(:, :, indxr),  xyVal] = ...
        blockRegisterMovie(data(:, :, indxr), ops.xyMask, dsall(indxr,:,:));
    
    ix0 = ix0 + Nbatch;
    xyValid = xyValid & xyVal;
end