function [icl, iSortFeat] = get_icl(icl, Nk, Ly, nFeat)

my = Inf * ones(Nk,1);
mx = Inf * ones(Nk,1);
icl = cell(Nk,1);

for i = 1:Nk
    ix = find(iclust==i);
    icl{i} = ix;
    my(i,1) = median(rem(ix-1, Ly)+1);
    mx(i,1) = median(ceil(ix/Ly));
end
ys = my * ones(1,Nk);
xs = mx * ones(1,Nk);

ds = (ys - ys').^2 + (xs - xs').^2;
[~, iSortFeat] = sort(ds, 1, 'ascend');
iSortFeat = iSortFeat(1:nFeat, :);
