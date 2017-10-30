function ds = RemoveBadShifts(ds0)

ds = ds0;

%%
ds_std = std(ds);
ds_mu  = mean(ds);
for j = 1:2
    inds = (ds(:,j) > ds_mu(j) + 5*ds_std(j)) | (ds(:,j) < ds_mu(j) - 5*ds_std(j));
    ds(inds,j) = NaN;
    if isnan(ds(1,j))
        ds(1,j) = ds(find(~isnan(ds(:,j)),1),j);
    end
    idx = ~isnan(ds(:,j));
    dn  = ds(idx,j);
    % take mean of point before and after nan
    dnm = (dn([1:end-1]) + dn([2:end]))/2;
    dnm = [dnm; dnm(end)];
    idx = cumsum(idx);
    idj = find(isnan(ds(:,j)));
    ds(idj,j) = dnm(idx(idj));
end
%%
clf
plot(ds0);
hold all;
plot(ds);
drawnow;
