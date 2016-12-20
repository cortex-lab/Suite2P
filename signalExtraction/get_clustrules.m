function clustrules = get_clustrules(clustrules)

clustrules.npix_fraclow             = getOr(clustrules, {'npix_fraclow'}, 1/2);
clustrules.npix_frachigh            = getOr(clustrules, {'npix_frachigh'}, 2);
clustrules.diameter                 = getOr(clustrules, {'diameter'}, 10);
clustrules.MinNpix                  = round(pi/4 * clustrules.diameter^2 *clustrules.npix_fraclow);
clustrules.MaxNpix                  = round(pi/4 * clustrules.diameter^2 *clustrules.npix_frachigh);
clustrules.Compact                  = getOr(clustrules, {'Compact'}, 1.1);
