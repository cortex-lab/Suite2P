% computes "footprint" of ROI and populates stat with cell locations
function stat = getFootprint(ops, codes, Ucell, mPix, mLam, stat)

d0   = ceil(ops.diameter); % expected cell diameter
dx = repmat([-2*d0:2*d0], 4*d0+1, 1);
dy = dx';

rs = (dx.^2 + dy.^2).^.5;
dx = dx(rs<=2*d0);
dy = dy(rs<=2*d0);
rs = rs(rs<=2*d0);

frac = [.5]; %[0.15 0.25 0.33 0.5 .75];
%%
Ny = length(ops.yrange);
Nx = length(ops.xrange);

% find maximum contamination distance for each ROI
for j = 1:size(mPix,2)
    ipos = find(mPix(:,j)>0 & mLam(:,j)>1e-10);
    ipix = mPix(ipos,j);
    
    stat(j).ipix    = ipix;
    [ypix, xpix]    = ind2sub([Ny Nx], ipix);
    stat(j).ypix    = ypix;
    stat(j).xpix    = xpix;
    stat(j).lam     = mLam(ipos,j);
    stat(j).lambda  = mLam(ipos,j);
    stat(j).npix    = numel(ipix);
    stat(j).med     = [median(stat(j).ypix) median(stat(j).xpix)];
    stat(j).neuropilCoefficient = 0.7;
    stat(j).baseline            = 0;
    
    % compute "footprint" length of ROI
    y0 = ceil(stat(j).med(1));
    x0 = ceil(stat(j).med(2));
    
    % pixels within FOV and within radius d0 of center
    ivalid = find((x0 + dx)>=1 & (x0 + dx)<=Nx & (y0 + dy)>=1 & (y0 + dy)<=Ny);
    ipix = (y0+dy(ivalid)) + (x0 + dx(ivalid)-1) * Ny;
    % weight of cell onto each pixel across components
    proj = codes(j,:) * Ucell(:, ipix); 
    % mean of radius of pixels with weights > frac
    stat(j).footprint = mean(rs(ivalid(proj>max(proj) * frac))) /ops.diameter;
    
end

mfoot = median([stat.footprint]);

for j = 1:size(mPix,2)
    stat(j).footprint = stat(j).footprint/mfoot;
    
%     stat(j).negFoot = max(0, 1 - stat(j).footprint/mfoot);
%     stat(j).posFoot = max(0, stat(j).footprint/mfoot -1);
end
