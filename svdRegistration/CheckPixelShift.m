function [cX,ix] = CheckPixelShift(A,B,isGPU)

if isGPU==1
    A = gpuArray(A);
    B = gpuArray(B);
end

% compute mean shift
[cc,mfactor] = regSVDcorr(A,B);

ccA = sum(cc.^2 .* ...
          repmat(permute(mfactor,[3 2 1]),size(cc,1),size(cc,2),1),3);

[cX,ixMax] =  max(ccA(:));

[ix(1),ix(2)] = ind2sub(size(ccA),ixMax);

if isGPU==1
    ix = gather(ix);
end