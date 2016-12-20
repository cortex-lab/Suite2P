function [caf, baselines, sn] = estimate_baseline(caCorrected, ops)

fs  = ops.fs;
trange = getOr(ops, 'trange', 60); % range of timepoints in seconds

caBase =  my_conv2(caCorrected, max(3, ceil(fs)), 1);
caBase =  my_min( caBase, ceil(fs*trange/2), 1);
caBase = -my_min(-caBase, ceil(fs*trange/2), 1);

caf = caCorrected - caBase;

[NT, NN] = size(caf);

sd = std(caf, 1, 1);
b0 = mean(caf,1);  
b1 = linspace(-1, 1, 20);

score = zeros(NN, length(b1));
for i = 1:length(b1)
    dF = bsxfun(@minus, caf, b0 + b1(i)*sd);
    dF = abs(dF);
    
    score(:,i) = mean(log(1+dF), 1);
end

[~, imax] = min(score, [], 2);

baselines = b0 + b1(imax).*sd;

caf  = bsxfun(@minus, caf, baselines);

baselines = caCorrected - caf;

fmin = .25;
[Pxx, fs] = pwelch(caf, [],[],[], 1);
sn = sqrt(exp(mean(log(Pxx(fs>fmin,:)/2), 1)));
caf = bsxfun(@rdivide, caf, sn);