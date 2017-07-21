% m1 = zstack yx x patches x depth (Z)
% m2 = yx     x patches
% iyxz = (Y X Z) best for each plane
function [iyxz,cbest]  = xcorr_patches(m1, m2)

cc        = real(ifft(ifft(bsxfun(@times, m1, conj((m2))),[],1),[],2));
czc       = squeeze(max(reshape(cc,[],size(cc,3),size(cc,4)),[],1));
[cbest,iz]    = max(czc, [], 2);
cxc       = max(reshape(cc,[],size(cc,3),size(cc,4)),[],3);
[~,ix]    = max(cxc, [], 1);
[iyc,ixc] = ind2sub([size(cc,1) size(cc,2)], ix);
ivec      = FindRegInds([iyc(:) ixc(:)], [size(cc,1) size(cc,2)]);

iz = squeeze(iz);
iz = iz - ceil(size(m1,4)/2);

iyxz = [ivec iz(:)];