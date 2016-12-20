function h = buildHue(h)

iclust          = h.dat.res.iclust;
H0              = zeros(h.dat.cl.Ly, h.dat.cl.Lx);
H0(iclust>0)    = h.dat.cl.rands(iclust(iclust>0));

h.dat.img1.H       = reshape(H0, h.dat.cl.Ly, h.dat.cl.Lx);
h.dat.img2.H       = h.dat.img1.H;
