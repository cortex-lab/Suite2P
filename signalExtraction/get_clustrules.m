function clustrules = get_clustrules(clustrules)

clustrules.npix_fraclow             = getOr(clustrules, {'npix_fraclow'}, 1/4);
clustrules.npix_frachigh            = getOr(clustrules, {'npix_frachigh'}, 20);
clustrules.diameter                 = getOr(clustrules, {'diameter'}, 10);
clustrules.MinNpix                  = round(pi/4 * clustrules.diameter^2 *clustrules.npix_fraclow);
clustrules.MaxNpix                  = round(pi/4 * clustrules.diameter^2 * clustrules.npix_frachigh);
clustrules.Compact                  = getOr(clustrules, {'Compact'}, 2);
clustrules.parent                   = getOr(clustrules, {'parent'}, []);
clustrules.parent.minPixRelVar      = getOr(clustrules.parent, {'minPixRelVar'}, 1/10);
clustrules.parent.PixelFractionThreshold     = getOr(clustrules.parent, {'PixelFractionThreshold'}, 0.5);
clustrules.parent.MaxRegions        = getOr(clustrules.parent, {'MaxRegions'}, 10);
