function h = splitROIleftright(h)
h.dat.cl.iscell = h.dat.cl.Mrs <h.dat.res.Mrs_thresh & ...
    h.dat.cl.npix <h.dat.cl.npix_high & h.dat.cl.npix >h.dat.cl.npix_low;

h.dat.cl.iscell = h.dat.cl.iscell  & (h.dat.cl.nreg      <h.dat.cl.nreg_max);
h.dat.cl.iscell = h.dat.cl.iscell  & (h.dat.cl.npix_par  <h.dat.cl.npix_par_max);
h.dat.cl.iscell = h.dat.cl.iscell  & (h.dat.cl.npix_res  <h.dat.cl.npix_res_max);
h.dat.cl.iscell = h.dat.cl.iscell  & (h.dat.cl.mrs_parent<h.dat.cl.mrs_parent_max);
h.dat.cl.iscell = h.dat.cl.iscell  & (h.dat.cl.VperPix   >h.dat.cl.VperPix_min);

h.dat.cl.iscell = double(h.dat.cl.iscell);
% overwrite manual selections
h.dat.cl.iscell(h.dat.cl.manual>1e-3) = 1;
h.dat.cl.iscell(h.dat.cl.manual<-1e-3) = 0;

h.dat.cl.k1 = reshape(h.dat.cl.iscell(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);

