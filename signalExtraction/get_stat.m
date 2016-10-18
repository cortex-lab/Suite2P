function stat = get_stat(res)

Ly = res.Ly;
Lx = res.Lx;

xs = repmat(1:Lx, Ly, 1);
ys = repmat((1:Ly)', 1, Lx);

xlx         = repmat(-ceil(Lx/2):1:ceil(Lx/2), 2*ceil(Lx/2)+1, 1);
rgrid       = sqrt(xlx.^2 + xlx'.^2);
rgridsort   = sort(rgrid(:), 'ascend');

Nk = numel(unique(res.iclust));
for k = 1:Nk
    ipix = find(res.iclust==k);    
    x0 = xs(ipix); y0 = ys(ipix);
    
    rs = ((x0 - median(x0)).^2 + (y0 - median(y0)).^2).^.5;        
    stat(k).mrs     = median(rs);
    stat(k).npix    = numel(ipix);
    stat(k).mrs0    = median(rgridsort(1:stat(k).npix));
    stat(k).med     = [median(y0) median(x0)];
    stat(k).ipix    = ipix;
    stat(k).lambda  = res.M(ipix);
    


%     stat(k).iscell = stat(k).mrs/stat(k).mrs0<1.3 & ...
%         stat(k).npix>50 & stat(k).npix<300;    
end