% registers frames using offsets dsall (rigid registration)
% loops over batches
function [dreg, xyValid] = nonrigidMovie(data, ops, dsall, xyValid)

[dreg,  xyVal] = ...
    nonrigidRegFrames(data, ops.xyMask, dsall);
    
   