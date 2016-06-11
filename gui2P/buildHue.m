function h = buildHue(h)
h.dat.img1.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
h.dat.img2.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
