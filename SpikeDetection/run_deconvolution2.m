% function [dcell, Ffr, kernel] = run_deconvolution2(Ff, Fneu, f0, kernel)
isroi = [stat.mrs]./[stat.mrs0]<1.2 & [stat.npix]>20 & [stat.npix]<200;
Ff = Fcell{1}(isroi, :);
Fneu = FcellNeu{1}(isroi, :);
load('D:\CODE\MariusBox\SpikeDetection\kernel.mat')
f0 = 3;

Ff   = Ff';
Fneu = -Fneu';

sd = std(Ff - my_conv2(Ff, 2, 1), [], 1);
Ff = 2 * Ff ./ repmat(1e-5 + sd, size(Ff,1), 1);
Fneu = 2 * Fneu ./ repmat(1e-5 + sd, size(Ff,1), 1);

Params = [1 3 3 3000]; %type of deconvolution, Th, Thi(nner loop), maxiter

kernel = interp1(1:numel(kernel), kernel, ...
    linspace(1, numel(kernel), ceil(f0/3 * numel(kernel))));
kernel = normc(kernel');
%%
clear dcell
tic
[NT, NN] = size(Ff);
Ffr = zeros(size(Ff));
for icell = 1:size(Ff,2)
    [dcell{icell}, Ffr(:, icell)]= deconv_single_cell(Ff(:, icell), ...
        Fneu(:,icell), Params, kernel);    
    plot(dcell{icell}.kernel)
    pause
end
toc
%%
dcell{1}.kernel = kernel;

Ffr = 1/2 * Ffr .* repmat(1e-5 + sd, size(Ff,1), 1);
%%
cneu = 1 + cellfun(@(x) x.B(3), dcell);

sk = skewness(my_conv2(Ff, 2, 1), [], 1);

plot(sk, cneu, 'o')