function stat = getFootprint(ops, codes, Ucell, mPix, mLam, mLam0, stat)

d0   = ceil(ops.diameter); % expected cell diameter
dx = repmat([-2*d0:2*d0], 4*d0+1, 1);
dy = dx';

rs = (dx.^2 + dy.^2).^.5;
dx = dx(rs<=2*d0);
dy = dy(rs<=2*d0);
rs = rs(rs<=2*d0);

frac = [.5]; %[0.15 0.25 0.33 0.5 .75];
%%
[nSVD, Ly, Lx] = size(Ucell);
% find maximum contamination distance for each ROI
for j = 1:size(mPix,2)
    ipos = find(mPix(:,j)>0 & mLam(:,j)>1e-10);
    ipix = mPix(ipos,j);
    
    stat(j).ipix    = ipix;
    stat(j).xpix    = ceil(ipix/Ly);
    stat(j).ypix    = rem(ipix-1, Ly)+1;
    stat(j).lam     = mLam(ipos,j);
    stat(j).lambda  = mLam0(ipos,j);
    stat(j).npix    = numel(ipix);
    stat(j).med     = [median(stat(j).ypix) median(stat(j).xpix)];
    
    % compute "footprint" length of ROI
    y0 = ceil(stat(j).med(1));
    x0 = ceil(stat(j).med(2));
    
    ivalid = find((x0 + dx)>=1 & (x0 + dx)<=Lx & (y0 + dy)>=1 & (y0 + dy)<=Ly);
    
    ipix = (y0+dy(ivalid)) + (x0 + dx(ivalid)-1) * Ly;
    proj = codes(j,:) * Ucell(:, ipix);
    
    stat(j).footprint = mean(rs(ivalid(proj>max(proj) * frac))) /ops.diameter;
    
end

mfoot = median([stat.footprint]);

for j = 1:size(mPix,2)
    stat(j).footprint = stat(j).footprint/mfoot;
    
%     stat(j).negFoot = max(0, 1 - stat(j).footprint/mfoot);
%     stat(j).posFoot = max(0, stat(j).footprint/mfoot -1);
end