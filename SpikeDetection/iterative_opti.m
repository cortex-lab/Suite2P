%% ITeRATIVE OPTIMIZATION OF PARAMETERS
addpath('D:\CODE\MariusBox\SpikeDetection')
%
find_kernel;
%%
icell = 10;

kernel = kernel/sum(kernel.^2).^.5 ;
F = dcell(icell).F;
% F = F - 2.5;

% F = my_conv(F', 1)';
sd = median(abs(F - median(F)));

%
Finp = (F - coefs(end))/sd;
Finp = double(Finp);
kernel = double(kernel);

% nt0 = numel(kernel);
npad = 250;

F1 = cat(1, zeros(npad,1), Finp, zeros(npad,1));

kerneli = kernel;
F1i = F1;
nt0 = numel(kernel);

taus = [1 2 4 8 16 32 64];
Nbasis = numel(taus);
kerns = zeros(nt0, Nbasis);
for i = 1:Nbasis
    kerns(:,i) = exp(-[1:nt0]/taus(i));
end
%%
NT = numel(F1i);
for iter = 1:3
    [F0, c, st] = deconvolveL0(F1i, kerneli, 'nonneg', 3);
    
    fs = zeros(numel(F1i), 1);
    fs(st+1) = c;
    
    X = zeros(NT, Nbasis+1);
    for i = 1:Nbasis
        X0 = conv(kerns(:,i), fs);
        X(:,i) = X0(1:NT);
    end
    
    X(:, Nbasis+1)= 1;
    
    if iter>1
        X(:, Nbasis+1) = [];
    end
    Cv = X' * X;
    B = Cv\(X' * F1);
    
    if iter==1
        F1i = F1 - B(Nbasis+1);
    end
    kerneli = kerns * B(1:Nbasis);
    kerneli = kerneli/sum(kerneli.^2).^.5 ;
    
    figure(2)
    plot(kerneli)
    drawnow
    numel(st)
    
    
    F0 = F0(npad+1:end-npad);
    st = st - npad;
    c(st<1 | st>numel(F)) = [];
    st(st<1 | st>numel(F)) = [];
    %
    hd = histc(dcell(icell).tspikes, 1:1:numel(F));
    fs = zeros(numel(F), 1);
    fs(st) = c;
    
    sig = [1 3 5 7 11 15 19 25 31 41 61 121];
    for i = 1:length(sig)        
        CC = corrcoef(smooth(fs(1:end-1)', sig(i) + 2), smooth(hd(2:end)', sig(i)));
        cfit(i,1) = CC(1,2).^2;
        
        CC = corrcoef(smooth(F(9:end)', sig(i) + 2), smooth(hd(1:end-8)', sig(i)));
        cfit(i,2) = CC(1,2).^2;        
    end
      
    K(:, iter) = kerneli;
    %     b(iter) = B(Nbasis + 1);
    %     disp(b(iter))
end
%%
s1 = fs;
s2 = hd;
s1 = zscore(s1, 1, 1)/numel(s1).^.5;
s2 = zscore(s2, 1, 1)/numel(s2).^.5;
cc = real(ifft(fft(s1, [], 1) .* fft(s2(end:-1:1), [], 1)));
plot(fftshift(cc))
%%
figure(1)
Frec = conv(fs, kerneli);
plot(Frec)
hold all
plot(Finp(2:end))

plot(-10 + 4 * fs/mean(fs(fs>1e-4)), 'Linewidth', 2, 'Linewidth', 1)
plot(-10  - 4*hd(2:end) , 'Linewidth', 2, 'Linewidth', 1)
hold off



