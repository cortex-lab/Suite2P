function [kernel, tau, Ball] = get1Dkernel(caHigh, fs)

Tlag = 1;

caHigh = double(caHigh);
% caHigh = bsxfun(@minus, caHigh, mean(caHigh,1));
% [b1, a1] = butter(3, 2 * 1/300 * 1/fs, 'high');
% caHigh   = filtfilt(b1, a1, caHigh);

[NT NN] = size(caHigh);

dataPred = caHigh(Tlag+1:NT,:);
err = zeros(NT, NN);

taus = exp(linspace(log(fs/10), log(2*fs), 20));

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
    Ball(:,j) = B;
        
    dE(j) = mean(err(:).^2);
end

[~, imin] = min(dE);
p = exp(-1/taus(imin));
%%

As = exp(-1./taus) + Ball(1,:);
tau = -1./log(As);

% disp(tau)

tau = tau(imin);

%%
kernel  = exp(-[0:1:ceil(5*tau)]'/tau);

% baselines = Ball(2,imin) / (1 - 1/(1-p) * Ball(1,imin));
