function mLam = getConnected(mLam, rs)

mLam1 = zeros(size(rs));
mLam1(rs<=0) = mLam>0;

rs0 = rs;
rs0(mLam1<1e-10) = Inf;
[~, imin] = min(rs0(:));

mask = zeros(size(rs));
mask(imin) = 1;

for i = 1:ceil(size(mask,1)/2)
    mask = my_max(mask, 1, 1) .* mLam1;
    mask = my_max(mask, 1, 2) .* mLam1;
end

mLam = mLam.*mask(rs<=0);