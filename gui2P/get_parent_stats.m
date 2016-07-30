function h = get_parent_stats(h)

h.dat.stat(1).V = [];
for i = 1:length(h.dat.stat)
    if ~isfield(h.dat.stat(i), 'V') || isempty(h.dat.stat(i).V)
        h.dat.stat(i).V         = sum([h.dat.stat(i).region.V]);
        h.dat.stat(i).mrs = min(h.dat.stat(i).mrs, 1e4);
       % this is a parent region 
        h.dat.stat(i).parent    = i;
        h.dat.stat(i).VperPix   = h.dat.stat(i).V/h.dat.stat(i).npix;
        h.dat.stat(i).npix_res  = numel(h.dat.stat(i).ipix);
        h.dat.stat(i).nregions  = numel(h.dat.stat(i).region);
    end
end

nreg = [h.dat.stat.nregions];
npix_res = [h.dat.stat.npix_res];
npix = [h.dat.stat.npix];

VperPix = [h.dat.stat.VperPix];
mrs = [h.dat.stat.mrs]./[h.dat.stat.mrs0];
iparent = [h.dat.stat.parent];

h.dat.cl.nreg       = nreg(iparent);
h.dat.cl.npix_res   = npix_res(iparent);
h.dat.cl.npix_par   = npix(iparent);
h.dat.cl.VperPix    = VperPix(iparent);
h.dat.cl.mrs_parent = mrs(iparent);
