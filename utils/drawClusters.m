function [r, iclust] = drawClusters(r, mPix, mLam, Ly, Lx, I)

icell = size(mPix,2);
iclust = (icell+1) * ones(Ly, Lx);

lam      = zeros(Ly, Lx);
for i = icell:-1:1
    ipos = find(mPix(:,i)>0);
    ipix = lam(mPix(ipos,i))+1e-4 < mLam(ipos,i);
    iclust(mPix(ipos(ipix),i)) = i;
    lam(mPix(ipos(ipix),i))      = mLam(ipos(ipix),i);
end

if nargin>5
    lam = I;
end
%
r = cat(1, r, rand(icell+1-numel(r), 1));
Sat = ones(Ly, Lx);

V = max(0, min(.75 * reshape(lam, Ly, Lx)/mean(lam(lam>1e-10)), 1));
H = reshape(r(iclust), Ly, Lx);
rgb_image = hsv2rgb(cat(3, H, Sat, V));
imagesc(rgb_image)
axis off