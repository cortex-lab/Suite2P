% find bad frames and ignore them during registration
function badi = getOutliers(ops)

Corr = ops.CorrFrame;
%%
ops.nSDbadregCorr   = getOr(ops, 'nSDbadregCorr', 5);
ops.nTbadregWindow  = getOr(ops, 'nTbadregWindow', 30);

ind = isnan(Corr);
x = 1:length(Corr);
vq = interp1(x(~ind),Corr(~ind),x(ind),'pchip');
Corr(ind) = vq;

mCorr = medfilt1(Corr, ops.nTbadregWindow);

lCorr = log(max(1e-6, Corr)) - log(max(1e-6, mCorr));

sd = std(lCorr);

badi = find(lCorr < -ops.nSDbadregCorr*sd);