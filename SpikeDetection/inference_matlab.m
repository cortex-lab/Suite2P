%% ITeRATIVE OPTIMIZATION OF PARAMETERS
addpath('D:\CODE\MariusBox\SpikeDetection')

icell = 1;
%
find_kernel;
%%
kernel = kernel/sum(kernel.^2).^.5 ;
F = dcell(icell).F;
F = F - 2.5;

sd = median(abs(F - median(F)));

%
Finp = (F - coefs(end))/sd;
Finp = double(Finp);
kernel = double(kernel);

% nt0 = numel(kernel);
npad = 250;

F0 = cat(1, zeros(npad,1), Finp, zeros(npad,1));
%
[F0, c, st] = deconvolveL0(F0, kernel, 'nonneg', 3);


F0 = F0(npad+1:end-npad);
st = st - npad;
c(st<1 | st>numel(F)) = [];
st(st<1 | st>numel(F)) = [];
%%
hd = histc(dcell(icell).tspikes, 1:1:numel(F));
fs = zeros(numel(F), 1);
fs(st) = c;

sig = 5;

CC = corrcoef(my_conv(fs(1:end-1)', sig), my_conv(hd(2:end)', sig));
cfit(1,1) = CC(1,2).^2;

CC = corrcoef(my_conv(F(9:end)', sig), my_conv(hd(1:end-8)', sig));
cfit(1,2) = CC(1,2).^2;
nanmean(cfit,1)

clf
Frec = conv(fs, kernel);
plot(Frec)
hold all
plot(Finp(2:end))


plot(-10 + 4 * fs/mean(fs(fs>1e-4)), 'Linewidth', 2, 'Linewidth', 1)
plot(-10  - 4*hd(2:end) , 'Linewidth', 2, 'Linewidth', 1)


%%
s1 = fs;
s2 = hd;
s1 = zscore(s1, 1, 1)/numel(s1).^.5;
s2 = zscore(s2, 1, 1)/numel(s2).^.5;
cc = real(ifft(fft(s1, [], 1) .* fft(s2(end:-1:1), [], 1)));
plot(fftshift(cc))

