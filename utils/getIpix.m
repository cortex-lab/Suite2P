function [ipix goodi] = getIpix(ind, dx, dy, Lx, Ly)

iy = rem(ind-1, Ly) + 1;
ix = ceil(ind/Ly);

yrange = iy + dy(:);
xrange = ix + dx(:);

badi = yrange<1 | yrange>Ly | xrange<1 | xrange>Lx;
goodi = find(~badi);
yrange = yrange(goodi);
xrange = xrange(goodi);

ipix = yrange + (xrange-1)*Ly;