function out = gpuBlockSmallXtY(X, Y)

g = gpuDevice;
bytesGPU = g.AvailableMemory;
infoX = whos('X');
infoY = whos('Y');
[hX, wX] = size(X);
[hY, wY] = size(Y); % hY should be equal to hX
% the dimensions of out will be wX-by-wY
bytesPerRowX = infoX.bytes/hX;
bytesPerRowY = infoY.bytes/hY;
% theoretically, to keep the result and {X', X, and Y} in GPU memory
% bytesRequired = bytesPerRowX*wY + batchSize*(2*bytesPerRowX+bytesPerRowY);
batchSize = floor((bytesGPU-bytesPerRowX*wY)/(2*bytesPerRowX+bytesPerRowY));
batchSize = floor(batchSize/8); % just to be on the safe side
% also, interestingly with smaller batches the code runs slightly faster
nBatches = ceil(hX/batchSize);

startIdx = 1:batchSize:hX;
endIdx = startIdx + batchSize-1;
endIdx = min(endIdx, hX);

out = gpuArray(zeros([wX, wY], 'like', X));
for iBatch = 1:nBatches
%     fprintf('Batch %d/%d\n', iBatch, nBatches);
    batchX = gpuArray(X(startIdx(iBatch):endIdx(iBatch), :));
    batchY = gpuArray(Y(startIdx(iBatch):endIdx(iBatch), :));
    out = out + batchX'*batchY;

% MK - this is an annoying line of code, but it helps Matlab to do memory
% management on the GPU properly. Otherwise I get OUT OF MEMORY errors
% maybe it is just my GPU is a bit unstable
    wait(gpuDevice);
end

out = gather(out);

return
