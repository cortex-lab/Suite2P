function [sp, dcell] = deconvolution_standalone2(ops, ca, neu, gt)
% takes as input the calcium and (optionally) neuropil traces,  
% both NT by NN (number of neurons).
% outputs a cell array dcell containing spike times (dcell.st) and amplitudes
% (dcell.c). dcell.B(3) is the neuropil contamination coefficient. 
% dcell.B(2) is an estimate of the baseline. 

% this version also estimates a single timescale

% specify in ops the following options, or leave empty for defaults
%       fs = sampling rate
%       recomputeKernel = whether to estimate kernel from data
%       sensorTau  = timescale of sensor, if recomputeKernel = 0
% additional options can be specified (mostly for linking with Suite2p, see below).

ops.imageRate    = getOr(ops, {'imageRate'}, 30); % total image rate (over all planes)
ops.nplanes      = getOr(ops, {'nplanes'}, 1); % how many planes at this total imaging rate
ops.sensorTau    = getOr(ops, {'sensorTau'}, 2); % approximate timescale in seconds
ops.sameKernel   = getOr(ops, {'sameKernel'}, 1); % 1 for same kernel per plane, 0 for individual kernels (does not work right now)
ops.maxNeurop    = getOr(ops, {'maxNeurop'}, Inf); % maximum allowed neuropil contamination coefficient. 
ops.recomputeKernel = getOr(ops, {'recomputeKernel'}, 1); % whether to estimate kernel from data
lam = getOr(ops, 'lam', 3);

% the kernel should depend on timescale of sensor and imaging rate
ops.fs           = getOr(ops, 'fs', ops.imageRate/ops.nplanes);
mtau             = ops.fs * ops.sensorTau; 

if nargin<3 || isempty(neu)
    neu = zeros(size(ca));
end

Params = [1 lam 1 2e4]; %parameters of deconvolution

% f0 = (mtau/2); % resample the initialization of the kernel to the right number of samples
kernel = exp(-[1:ceil(5*mtau)]'/mtau) ;
%
npad        = 250;
[NT, NN]    = size(ca);
coefNeu     = .8 * ones(1,NN); % initialize neuropil subtraction coef with 0.8

caCorrected = ca - bsxfun(@times, neu, coefNeu);

if ops.recomputeKernel
    tlag                     = 1;
    [kernel, mtau]           = estimateKernel(ops, caCorrected, tlag);
    
%     [kernel, mtau, coefNeu]  = estimateKernel(ops, ca - coefNeu * neu, tlag);
    
    fprintf('Timescale determined is %4.4f samples \n', mtau);
end
kernel = normc(kernel(:));
%%
kernelS     = repmat(kernel, 1, NN);
dcell       = cell(NN,1);

% running minimum subtraction here???


Fsort       = my_conv2(caCorrected, ceil(ops.fs), 1);


Fsort       = sort(Fsort, 1, 'ascend');
baselines   = Fsort(ceil(NT/20), :);

% determine and subtract the neuropil
F1 = caCorrected - bsxfun(@times, ones(NT,1), baselines);

% normalize signal
% sd   = 1/2 * std(F1 - my_conv2(F1, max(2, ops.fs/4), 1), [], 1);

sd   = 1/2 * 1/sqrt(2) * std(F1(2:end, :) - F1(1:end-1, :), [], 1);

F1   = bsxfun(@rdivide, F1 , 1e-12 + sd);


if 1
    % get new kernel
    kerns = exp(-bsxfun(@rdivide, [1:ceil(5*mtau)]',  mtau * [.05 .125 .25 .5 1 2]));
    spfilt = ones([size(kerns,2) size(F1)]);
    for j = 1:size(kerns,2)
        spfilt(j, :,:) = filter(kerns(:,j), 1, gt);
    end
    
    sts = zeros(size(spfilt,1));
    stF = zeros(size(spfilt,1), 1);
    for k = 1:size(spfilt,3)
        sts = sts + spfilt(:,:,k) * spfilt(:,:,k)';
        stF = stF + spfilt(:,:,k) * F1(:,k);
    end
    coefs = sts \ stF;
    kernel = normc(kerns * coefs);
    kernelS = repmat(kernel, 1, NN);

    plot(kernel)
    drawnow
end

sp = zeros(size(F1));

% run the deconvolution to get fs etc\
parfor icell = 1:size(ca,2)
    [sp(:,icell),dcell{icell}] = ...
        single_step_single_cell(F1(:,icell), Params, kernelS(:,icell), NT, npad,dcell{icell});
end

% rescale baseline contribution
for icell = 1:size(ca,2)
    dcell{icell}.c                      = dcell{icell}.c * sd(icell);
    dcell{icell}.baseline               = baselines(icell);
    dcell{icell}.neuropil_coefficient   = coefNeu(icell);
end

%%

