function h = splitROIleftright(h)
h.dat.cl.iscell = h.dat.cl.Mrs <h.dat.res.Mrs_thresh & ...
    h.dat.cl.npix <h.dat.cl.npix_high & h.dat.cl.npix >h.dat.cl.npix_low;

h.dat.cl.iscell = double(h.dat.cl.iscell);

% overwrite manual selections
h.dat.cl.iscell(h.dat.cl.manual>1e-3) = 1;
h.dat.cl.iscell(h.dat.cl.manual<-1e-3) = 0;
%%
IScell = zeros(h.dat.cl.Ly, h.dat.cl.Lx);
ix = h.dat.res.iclust>0;
IScell(ix) = h.dat.cl.iscell(h.dat.res.iclust(ix));
h.dat.cl.k1 = reshape(IScell, h.dat.cl.Ly, h.dat.cl.Lx);

