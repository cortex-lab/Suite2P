function [footPrint, nPix]= get_footprint(xs, Ly, Lx, ops)


[Mmax, imax] = max(xs, [], 2);

yp = rem(imax-1, Ly) + 1;
xp = ceil(imax/Ly);

Nk = size(xs,1);
footPrint = zeros(Nk, 1);
nPix = zeros(Nk, 1);

d = ops.diameter;
XS = repmat([-d*2:d*2], 4*d+1, 1);
YS = XS';

ds = (XS.^2 + YS.^2).^.5;

for j = 1:Nk
    yp1 = yp(j) + YS;
    xp1 = xp(j) + XS;
    
    badi = (yp1<1) | (yp1>Ly) | (xp1<1) | (xp1>Lx);
    yp1 = yp1(~badi);
    xp1 = xp1(~badi);
    allds = ds(~badi);
    
    ind = xs(j, yp1 + (xp1-1) * Ly) > Mmax(j)/2;
    
    footPrint(j) = mean(allds(ind));
    nPix(j) = numel(allds);
end

end