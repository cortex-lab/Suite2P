function [kernel, tau] = get1Dkernel(caHigh, Mtau0, Tlag, fs)

% caHigh = ca;
[b1, a1] = butter(3, 2 * 1/50 * 1/fs, 'high');
caHigh   = filtfilt(b1, a1, double(caHigh));

[NT NN] = size(caHigh);

dataPred = caHigh(Tlag+1:NT,:);
err = zeros(NT, NN);

taus = exp(linspace(log(1), log(3*fs), 10));

dE = zeros(length(taus), 1);
Ball = zeros(2, length(taus));
for j = 1:length(taus)
    Mtau = taus(j);
    
    kernel  = exp(-[0:1:ceil(5*Mtau)]'/Mtau);
    cafilt  = filter(kernel, 1, caHigh);
    cafilt  = cafilt(1:NT-Tlag,:);
    
    Y       = dataPred(:)/NT;
    X       = [cafilt(:) ones(numel(Y), 1)]/NT;
    B       = (X'*X)\(X'*Y);
    
    Ypred = X*B;
    
    err(1:NT-Tlag, :)   = reshape(Ypred-Y, size(dataPred));
    dErr                = filter(kernel, 1, flipud(B(1) * err));
    
    Ball(:,j) = B;
    dp(j) = mean(mean(dErr(2:NT-Tlag, :) .* cafilt(1:NT-Tlag-1, :)));
    
    
    dE(j) = mean(err(:).^2);
end

[~, imin] = min(dE);
p = exp(-1/taus(imin));
%%
tau2 = exp(linspace(log(0.5*fs), log(5*fs), 30));
As = exp(-1./tau2);

tau = fs * ones(1, numel(taus));
for i  = 1:numel(taus)
    ibest = find((As - exp(-1/taus(i))) .* As.^Tlag > Ball(1,i), 1);
    if ~isempty(ibest)
        tau(1,i) = -1./log(As(ibest));
    end
end
disp(tau)

% tau = tau(min(imin, numel(tau)));
tau = tau(imin);

%%
kernel  = exp(-[0:1:ceil(5*tau)]'/tau);