% computes compactness of ROIs and fits ellipses
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

%% neighbors for obtaining exterior of ROIs
Ly     = length(ops.yrange);
Lx     = length(ops.xrange);
ipts   = [1:Ly*Lx]';
ineigh = [ipts-1 ipts+1 ipts-Ly ipts+Ly];
ineigh(ineigh<1 | ineigh > Ly*Lx) = NaN;


%%

for j = 1:size(mPix,2)
    
    lam  = mLam(:,j);
        
    gpix = lam>1e-3;

    dd = d2p(gpix, gpix);

    
    stat(j).mrs(1) = mean(dd(:))/d0;     
    dd = d2p0(1:sum(gpix), 1:sum(gpix));
    stat(j).mrs0(1) = mean(dd(:))/di;    
    
    stat(j).cmpct = stat(j).mrs(1)/stat(j).mrs0(1);

   
    params = FitMVGaus(stat(j).ypix,stat(j).xpix, stat(j).lam);
    
    % save ellipse information
    stat(j).aspect_ratio = sqrt(max(params.eval) / min(params.eval));
    stat(j).ellipse      = params.xy;
    
end



%%