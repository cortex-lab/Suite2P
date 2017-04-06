% compute X*Y on GPU in batches to avoid memory errors
function out = gpuBlockXY(X, Y)

% this code is for a very tall X, and small Y

g = gpuDevice;
bytesGPU = g.AvailableMemory;
infoX = whos('X');
infoY = whos('Y');
[hX, wX] = size(X);
[hY, wY] = size(Y); % hY should be equal to wX
% the dimensions of out will be hX-by-wY
bytesPerRowX = infoX.bytes/hX;
bytesPerRowY = infoY.bytes/hY;
% theoretically, to keep the tmp result, X, and Y in GPU memory
% bytesRequired = bytesPerRowY*batchSize + bytesPerRowX*batchSize + infoY.bytes;
batchSize = floor((bytesGPU-infoY.bytes)/(bytesPerRowX+bytesPerRowY));
batchSize = floor(batchSize/8); % just to be on the safe side
% also, interestingly with smaller batches the code runs slightly faster
nBatches = ceil(hX/batchSize);

startIdx = 1:batchSize:hX;
endIdx = startIdx + batchSize-1;
endIdx = min(endIdx, hX);

out = zeros([hX, wY], 'like', X);
gpuY = gpuArray(Y);
for iBatch = 1:nBatches
%     fprintf('Batch %d/%d\n', iBatch, nBatches);
    batchX = gpuArray(X(startIdx(iBatch):endIdx(iBatch), :));
    tmp = batchX*gpuY;

% MK - this is an annoying line of code, but it helps Matlab to do memory
% management on the GPU properly. Otherwise I get OUT OF MEMORY errors
% maybe it is just my GPU is a bit unstable
    wait(g);
    out(startIdx(iBatch):endIdx(iBatch), :) = gather(tmp);
end

return
