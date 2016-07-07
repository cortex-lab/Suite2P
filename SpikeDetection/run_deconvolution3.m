function [dcell, Ffr] = run_deconvolution3(ops, Ff, Fneu, kernel)
% outputs a cell array of deconvolved spike times and amplitudes.
% Optionally output this in matrix form Ffr (very sparse). 

% the basis functions should depend on timescale of sensor and imaging rate
mtau             = ops.imageRate * ops.sensorTau/ops.nplanes; 

coefNeu = 0.8; % initialize neuropil subtraction coef with 0.8

sd   = std(Ff - my_conv2(Ff, 2, 1), [], 1);
Ff   = 2 * Ff ./ repmat(1e-5 + sd, size(Ff,1), 1);
Fneu = 2 * Fneu ./ repmat(1e-5 + sd, size(Ff,1), 1);

Params = [1 3 3 2e4]; %type of deconvolution, Th, Thi(nner loop), max Nspikes

f0 = (mtau/2); % resample the initialization of the kernel to the right number of samples
kernel = interp1(1:numel(kernel), kernel, ...
    linspace(1, numel(kernel), ceil(f0/3 * numel(kernel))));
kernel = normc(kernel');

npad = 250;
[NT, NN] = size(Ff);
F1 = cat(1, zeros(npad,NN), double(Ff - coefNeu * Fneu), zeros(npad,NN));
F1 = bsxfun(@minus, F1 , median(F1,1));

nt0 = numel(kernel);

taus = mtau * [1/4 1/2 1 2 4];
Nbasis = numel(taus);
kerns = zeros(nt0, Nbasis);
for i = 1:Nbasis
    kerns(:,i) = exp(-[1:nt0]/taus(i));
end
%

kernelS = repmat(kernel, 1, NN);

for iter = 1:10
    fprintf('iter %d... ', iter)
    parfor icell = 1:size(Ff,2)
        [kernelS(:, icell), F1(:,icell)] =...
            single_step_single_cell(ops, Ff(:,icell), F1(:, icell), Fneu(:,icell), Params, ...
            kernelS(:,icell), kerns, NT, npad);
    end
    if ops.sameKernel
        kernelS = repmat(median(kernelS,2), 1, NN);
    end
    
end
fprintf('finished, final extraction step... \n')
%%
dcell = cell(NN,1);
Ffr = zeros(size(Ff));
parfor icell = 1:size(Ff,2)
    [~, ~, dcell{icell}, Ffr(:, icell)] = single_step_single_cell(ops, ...
        Ff(:,icell), F1(:, icell), Fneu(:,icell), Params, ...
        kernelS(:,icell), kerns, NT, npad, dcell{icell});
end

% rescale baseline contribution
for icell = 1:size(Ff,2)
   dcell{icell}.B(2) = dcell{icell}.B(2) * sd(icell);
end

% rescale the deconvolved trace
Ffr = 1/2 * Ffr .* repmat(1e-5 + sd, size(Ff,1), 1);
%%
