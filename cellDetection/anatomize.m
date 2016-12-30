function stat = anatomize(ops, mPix, mLam, stat)

di = ops.diameter;
d0   = ceil(ops.diameter); % expected cell diameter

% data Ucell is nMaps by Ly by Lx

dx = repmat([-d0:d0], 2*d0+1, 1);
dy = dx';

rs = dx.^2 + dy.^2;
dx = dx(rs<=d0^2);
dy = dy(rs<=d0^2);

d2p = (bsxfun(@minus, dx, dx').^2 + bsxfun(@minus, dy, dy').^2).^.5;

xlx         = repmat(-ceil(2*d0):1:ceil(2*d0), 2*ceil(2*d0)+1, 1);
rgrid       = sqrt(xlx.^2 + xlx'.^2);
[rgridsort, isort]  = sort(rgrid(:), 'ascend');
xlxt        = xlx';

d2p0 = (bsxfun(@minus, xlx(:), xlx(:)').^2 + bsxfun(@minus, xlxt(:), xlxt(:)').^2).^.5;
d2p0 = d2p0(isort, isort);

%%
rd = zeros(size(mPix,2), 1);
rd0 = zeros(size(mPix,2), 1);


for j = 1:size(mPix,2)
    
    lam  = mLam(:,j);
        
    gpix = lam>1e-3;

    dd = d2p(gpix, gpix);
    stat(j).mrs(1) = mean(dd(:))/d0; %mean(ds(gpix));
    
%     ds = ((dx(gpix) - median(dx(gpix))).^2 + ...
%         (dy(gpix) - median(dy(gpix))).^2).^.5;
%     stat(j).mrs(2) = mean(ds(:))/di; %mean(ds(gpix));
    
    dd = d2p0(1:sum(gpix), 1:sum(gpix));
    stat(j).mrs0(1) = mean(dd(:))/di;    
    
    stat(j).cmpct = stat(j).mrs(1)/stat(j).mrs0(1);
%     stat(j).mrs0(2) = mean(rgridsort(1:sum(gpix)))/di;    
end



%%