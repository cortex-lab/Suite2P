function badi = getOutliers(ops)

Corr = ops.CorrFrame;
%%
ops.nSDbadregCorr   = getOr(ops, 'nSDbadregCorr', 5);
ops.nTbadregWindow  = getOr(ops, 'nTbadregWindow', 30);

mCorr = medfilt1(Corr, ops.nTbadregWindow);

sd = std(Corr - mCorr);

badi = find(Corr < mCorr - ops.nSDbadregCorr*sd);