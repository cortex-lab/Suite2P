% cut through z-stack using transformation matrix A
% Mimg = zstack
% Ly, Lx size of imaged plane
% zspread is number of planes above and below imaged plane to keep in
% z-stack
function MimgR = cutZstack(Mimg, A, Ly, Lx, zspread)

[Ny Nx Nz] = size(Mimg);

[Y X]      = ndgrid([1:Ly], [1:Lx]);

yxz        = round(A * [Y(:)'; X(:)'; ones(1, numel(Y))]);
yxz        = yxz';
yxz        = min(yxz, repmat(size(Mimg), size(yxz,1), 1));
yxz        = max(1, yxz);

% make a z-stack around aligned plane
if nargin > 4
    minZ   = min(yxz(:,3));
    maxZ   = max(yxz(:,3));
    MimgR  = NaN*zeros(Ly, Lx, 2*zspread + 1);
    zinds  = [-1*zspread + max(0, zspread-minZ + 1) : zspread - max(0, (maxZ+zspread)-Nz)];
    for j = [1:length(zinds)]
        zval = zinds(j) + yxz(:,3);
        ipix = sub2ind([Ny Nx Nz], yxz(:,1), yxz(:,2), zval);
        MimgR(:,:, j + max(0,zspread-minZ+1)) = reshape(Mimg(ipix), Ly, Lx);
    end
else
    ipix       = sub2ind([Ny Nx Nz], yxz(:,1), yxz(:,2), yxz(:,3));
    MimgR      = reshape(Mimg(ipix), Ly, Lx);
end
   