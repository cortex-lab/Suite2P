function [dcell, isroi] = run_deconvolution3(ops, dat, kernel)
% outputs a cell array of deconvolved spike times and amplitudes.
% Optionally output this in matrix form Ffr (very sparse).

% load the initialization of the kernel    
load(fullfile(ops.toolbox_path, 'SpikeDetection\kernel.mat'));

if isfield(dat, 'cl') && isfield(dat.cl, 'iscell')
    isroi = dat.cl.iscell;
else
    isroi = [dat.stat.mrs]./[dat.stat.mrs0]<dat.clustrules.Compact & ...
        [dat.stat.npix]>dat.clustrules.MinNpix & [dat.stat.npix]<dat.clustrules.MaxNpix;
end
isroi = logical(isroi);

% construct Ff and Fneu
Ff = [];
Fneu = [];
if isfield(dat, 'F')
    for j = 1:numel(dat.F.Fcell)
        Ff   = cat(1, Ff, dat.F.Fcell{j}(isroi, :)');
        Fneu = cat(1,Fneu, dat.F.FcellNeu{j}(isroi, :)');
    end
    flag = mean(sign(dat.F.FcellNeu{1}(:)))<0;
else
    for j = 1:numel(dat.Fcell)
        Ff   = cat(1, Ff, dat.Fcell{j}(isroi, :)');
        Fneu = cat(1,Fneu, dat.FcellNeu{j}(isroi, :)');
    end
    flag = mean(sign(dat.FcellNeu{1}(:)))<0;
end

if flag
    % then it means this was processed with old "model" option
    Fneu = -Fneu;
    Ff   = Ff + Fneu;
end

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
%
npad = 250;
[NT, NN] = size(Ff);
%%
B = zeros(3, NN);
B(2,:) = median(Ff - coefNeu * Fneu,1);
B(3,:) = coefNeu;

nt0 = numel(kernel);

taus = mtau * [1/8 1/4 1/2 1 2 4];
Nbasis = numel(taus);
kerns = zeros(nt0, Nbasis);
for i = 1:Nbasis
    kerns(:,i) = exp(-[1:nt0]/taus(i));
end
%

kernelS = repmat(kernel, 1, NN);

maxNeurop = ops.maxNeurop;
tic
err = zeros(NN, 1);
for iter = 1:10
    fprintf('%2.2f sec, iter %d... ', toc, iter)
    parfor icell = 1:100 %size(Ff,2)
        [kernelS(:, icell), B(:,icell), err(icell)] = ...
            single_step_single_cell(maxNeurop, Ff(:,icell), B(:, icell), Fneu(:,icell), Params, ...
            kernelS(:,icell), kerns, NT, npad);
        
    end
%     mean(err(:))
    
    if ops.sameKernel
        kernelS = repmat(normc(nanmedian(kernelS(:, 1:100),2)), 1, NN);
    end
    kernelS = max(0, kernelS); % make sure kernel is positive
    kernelS  = normc(kernelS) ; % renormalize the kernel

    plot(kernelS(:,1:100))
    drawnow
end

fprintf('finished, final extraction step... \n')
%%
dcell = cell(NN,1);
for icell = 1:size(Ff,2)
    [~, ~, dcell{icell}] = ...
        single_step_single_cell(maxNeurop, Ff(:,icell), B(:,icell), Fneu(:,icell), Params, ...
        kernelS(:,icell), kerns, NT, npad, dcell{icell});
end

% rescale baseline contribution
for icell = 1:size(Ff,2)
   dcell{icell}.B(2) = dcell{icell}.B(2) * sd(icell);
end

%%

