function [cx,ix,ccA] = SVDPixelShift(A,B,nLocalPix)

numBlocks = size(A,3);
nX = size(A,1); npix = nX;
nY = size(A,2);

ccA = zeros(nLocalPix,nLocalPix,'single');

% compute correlation matrix
[cc,mfactor] = regSVDcorr(A,B);

cc0 = fftshift(fftshift(cc,1),2);
fftCent = [floor(nX/2)+1 floor(nY/2)+1];
dX = [-nLocalPix:nLocalPix]; dY = dX;
cc     = cc0(fftCent(1)+dX,fftCent(2)+dY,:);
cc     = ifftshift(ifftshift(cc,1),2);

% combine 10x10
ccA = sum(cc.^2 .* ...
          repmat(permute(mfactor,[3 2 1]),size(cc,1),size(cc,2),1),3);

% find max of correlation matrix
[cx,ixmax] = max(ccA(:));
[ix1,ix2] = ind2sub(size(ccA),ixmax);
ix = [ix1 ix2];

% shift ixs
ix = FindRegInds(ix,size(cc,1),size(cc,2));

