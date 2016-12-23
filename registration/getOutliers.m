function badi = getOutliers(ops)

Corr = ops.CorrFrame;
%%
ops.nSDbadregCorr   = getOr(ops, 'nSDbadregCorr', 5);
ops.nTbadregWindow  = getOr(ops, 'nTbadregWindow', 30);

mCorr = medfilt1(Corr, ops.nTbadregWindow);

lCorr = log(max(1e-6, Corr)) - log(max(1e-6, mCorr));

sd = std(lCorr);

badi = find(lCorr < -ops.nSDbadregCorr*sd);