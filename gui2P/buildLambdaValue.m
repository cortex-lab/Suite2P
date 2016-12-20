function h = buildLambdaValue(h)
h.dat.img1.V       = h.dat.img0.V .* h.dat.cl.k1;
h.dat.img2.V       = h.dat.img0.V .* (~h.dat.cl.k1);

iselect = h.dat.res.iclust==h.dat.F.ichosen;
h.dat.img1.V(iselect) = h.dat.cl.k1(iselect);
h.dat.img2.V(iselect) = ~h.dat.cl.k1(iselect);
