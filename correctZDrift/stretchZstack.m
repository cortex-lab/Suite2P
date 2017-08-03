% stretch z-stack so that it is at same angle in brain as multi-plane
% imaged planes
% zpos is position of each plane in zstack
% angle is in radians
function MimgR = stretchZstack(Mimg, ang)

[Ly Lx Lz] = size(Mimg);

dZ         = Ly * tan(ang);
[X Z]      = ndgrid([1:Lx], [1:Lz]);
X2         = X;

MimgR = zeros(size(Mimg),'single');

for i = 1:Ly
    Z2 = Z - (dZ/Ly)*(i-1);
    M0 = squeeze(Mimg(i,:,:));
    MimgR(i,:,:) = interp2([1:Lz],[1:Lx],M0,Z2,X2,'linear');
end

MimgR(isnan(MimgR(:))) = 0;