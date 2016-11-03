function out = gpuBlockSmallXtY(X, Y)

% X is small, Y is long (shallow and wide)

g = gpuDevice;
bytesGPU = g.AvailableMemory;
infoX = whos('X');
infoY = whos('Y');
[hX, wX] = size(X);
[hY, wY] = size(Y); % hX should be equal to hY
% the dimensions of 'out' will be wX-by-wY
bytesPerColumnY = infoY.bytes/wY;
% theoretically, to keep the tmp, X' and batchY in GPU memory
% bytesRequired = wX*batchSize + infoX.bytes+batchSize*bytesPerColumnY;
batchSize = floor((bytesGPU-infoX.bytes)/(bytesPerColumnY+wX));
batchSize = floor(batchSize/8); % just to be on the safe side
% also, interestingly with smaller batches the code runs slightly faster
nBatches = ceil(wY/batchSize);

startIdx = 1:batchSize:wY;
endIdx = startIdx + batchSize-1;
endIdx = min(endIdx, wY);

out = zeros([wX, wY], 'like', X);
out(end) = 1;
gpuX = gpuArray(X);
gpuXT = gpuX';
for iBatch = 1:nBatches
%     fprintf('Batch %d/%d\n', iBatch, nBatches);
    batchY = gpuArray(Y(:, startIdx(iBatch):endIdx(iBatch)));
    tmp = gpuXT*batchY;

% MK - this is an annoying line of code, but it helps Matlab to do memory
% management on the GPU properly. Otherwise I get OUT OF MEMORY errors
% maybe it is just my GPU is a bit unstable
    wait(g);
    out(:, startIdx(iBatch):endIdx(iBatch)) = gather(tmp);
end


return
