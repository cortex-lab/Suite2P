% compute covariance matrix of A on GPU in batches to avoid memory errors
function out = gpuBlockXtX(A)

g = gpuDevice;
bytesGPU = g.AvailableMemory;
infoA = whos('A');
[h, w] = size(A);
bytesPerRow = infoA.bytes/h;
% theoretically, to keep the result and two (A dn A') copies of data
% bytesRequired = bytesPerRow*w + 2*batchSize*bytesPerRow;
batchSize = floor((bytesGPU/bytesPerRow-w)/2);
batchSize = floor(batchSize/8); % just to be on the safe side
% also, interestingly with smaller batches the code runs slightly faster
nBatches = ceil(h/batchSize);

startIdx = 1:batchSize:h;
endIdx = startIdx + batchSize-1;
endIdx = min(endIdx, h);

out = gpuArray(zeros(w, 'like', A));
for iBatch = 1:nBatches
%     fprintf('Batch %d/%d\n', iBatch, nBatches);
    batchA = gpuArray(A(startIdx(iBatch):endIdx(iBatch), :));
    out = out + batchA'*batchA;
    
    % trying to debug OUT OF MEMORY issues on GPU
%     batchAT = batchA';
%     out = out + batchAT*batchA;
%     clear batchA batchAT;

% MK - this is an annoying line of code, but it helps Matlab to do memory
% management on the GPU properly. Otherwise I get OUT OF MEMORY errors
% maybe it is just my GPU is a bit unstable
    wait(gpuDevice);
end

out = gather(out);

return

%% this is a non-GPU version
batchSize = 2^16;
[h, w] = size(A);
out = zeros(w, 'like', A);
nBatches = ceil(h/batchSize);
startIdx = 1:batchSize:h;
endIdx = startIdx + batchSize-1;
endIdx = min(endIdx, h);

for iBatch = 1:nBatches
    fprintf('Batch %d/%d\n', iBatch, nBatches);
    batchA = A(startIdx(iBatch):endIdx(iBatch), :);
    out = batchA'*batchA + out;
end

