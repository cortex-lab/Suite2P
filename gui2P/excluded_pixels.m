function h = excluded_pixels(h)
h.dat.cl.excluded_pixels = h.dat.res.M<h.dat.cl.pixthresh_var;
h.dat.cl.excluded_pixels = reshape(h.dat.cl.excluded_pixels, h.dat.cl.Ly, h.dat.cl.Lx);