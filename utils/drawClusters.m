function [r, iclust, lam] = drawClusters(ops, r, mPix, mLam, Ly, Lx)

icell = size(mPix,2);
iclust = zeros(Ly, Lx);

lam      = zeros(Ly, Lx);
for i = icell:-1:1
    ipos = find(mPix(:,i)>0);
    ipix = lam(mPix(ipos,i))+1e-4 < mLam(ipos,i);
    iclust(mPix(ipos(ipix),i)) = i;
    lam(mPix(ipos(ipix),i))      = mLam(ipos(ipix),i);
end

% if nargin>5
%     lam = I;
% end

r = cat(1, r, rand(icell+1-numel(r), 1));

if ops.fig
%     figure(2)
    clf
    Sat = ones(Ly, Lx);
    H = zeros(Ly, Lx);
    
    H(iclust>0) = r(iclust(iclust>0));
    
    V = max(0, min(.75 * reshape(lam, Ly, Lx)/mean(lam(lam>1e-10)), 1));
    H = reshape(H, Ly, Lx);
    rgb_image = hsv2rgb(cat(3, H, Sat, V));
    imagesc(rgb_image)
    axis off
end

drawnow