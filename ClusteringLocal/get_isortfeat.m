function iSortFeat   = get_isortfeat(icl, Nk, Ly, nFeat)

ys = cellfun(@(x) mean(rem(x-1, Ly)+1), icl);
xs = cellfun(@(x) mean(ceil(x/Ly)), icl);

ys(isnan(ys)) = Inf;
xs(isnan(xs)) = Inf;

ds = bsxfun(@minus, ys, ys').^2 + bsxfun(@minus, xs, xs').^2;

[~, iSortFeat] = sort(ds, 1, 'ascend');

iSortFeat = iSortFeat(1:nFeat, :);
