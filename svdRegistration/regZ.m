function [cmax,ix,cc] = regZ(M1,M2)

[nX nY] = size(M1);
eps0 = single(1e-20);

m1 = fftn(M1);
m1 = m1./(abs(m1)+eps0);

m2 = fftn(M2);
m2 = m2./(abs(m2)+eps0);

cc = real(ifftn(m1 .* conj(m2)));

cx = []; ix = [];

[cmax,ixmax] = max(cc(:));
[ix1,ix2] = ind2sub(size(cc),ixmax);
  
ix = [ix1 ix2];

% shift ixs
nA = [nX nY];
ixind = ix > nA/2;

ix(ixind) = ix(ixind) - (nA(ixind)+1);
ix(~ixind) = ix(~ixind)-1;


% lCorr = 50;
% [yClipRef, xClipRef] = ndgrid(-2:2, -2:2);
% xClipRef = xClipRef(:);
% yClipRef = yClipRef(:);
% nClipPixels = numel(xClipRef);
% xCorrRef = [(usFac*lx - lCorr + 1):usFac*lx 1:(lCorr + 1)];
% yCorrRef = [(usFac*ly - lCorr + 1):usFac*ly 1:(lCorr + 1)];
% % weighted average for subpixel shifts
% iy = min(max(iy, 3), 2*lCorr - 1);
%     ix = min(max(ix, 3), 2*lCorr - 1);
%     clipX = bsxfun(@plus, xClipRef', ix);
%     clipY = bsxfun(@plus, yClipRef, iy);
%     clipF = reshape(repmat(1:size(clipX, 3), nClipPixels, 1), [], 1);
%     cczoom = reshape(...
%       gather(corrClip(sub2ind(size(corrClip), clipY(:), clipX(:), clipF))),...
%       nClipPixels, 1, []);
%     bcorr = sum(cczoom, 1);
%     cczoom = bsxfun(@rdivide, cczoom, bcorr);
%     ix = ix + sum(bsxfun(@times, xClipRef, cczoom), 1);
%     iy = iy + sum(bsxfun(@times, yClipRef, cczoom), 1);
end