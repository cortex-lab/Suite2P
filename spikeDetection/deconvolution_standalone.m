function [sp, dcell] = deconvolution_standalone(ops, ca, neu)
% takes as input the calcium and (optionally) neuropil traces,  
% both NT by NN (number of neurons).

% outputs an array of deconvolved spike times and amplitudes
% also outputs a cell array dcell containing the same information
% (dcell.c). dcell.B(3) is the neuropil contamination coefficient. 
% dcell.B(2) is an estimate of the baseline. 

ops.imageRate    = getOr(ops, {'imageRate'}, 30); % total sampling rate (over all planes)
ops.nplanes      = getOr(ops, {'nplanes'}, 1); % how many planes at this total imaging rate
ops.sensorTau    = getOr(ops, {'sensorTau'}, 2); % approximate timescale in seconds
ops.sameKernel   = getOr(ops, {'sameKernel'}, 1); % 1 for same kernel per plane, 0 for individual kernels (does not work right now)
ops.maxNeurop    = getOr(ops, {'maxNeurop'}, Inf); % maximum allowed neuropil contamination coefficient. 
ops.recomputeKernel = getOr(ops, {'recomputeKernel'}, 1); % whether to estimate kernel from data

% the sampling rate can also be specified directly
ops.fs           = getOr(ops, 'fs', ops.imageRate/ops.nplanes);
% initialize timescale of kernel
mtau             = ops.fs * ops.sensorTau; 
kernel = exp([-1:ceil(5*mtau)]'/mtau);

% fill up neu with zeros if absent
if nargin<3 || isempty(neu)
    neu = zeros(size(ca));
end

% size of ca array
[NT, NN]    = size(ca);

% initialize neuropil subtraction coef with 0.8
coefNeu = .8 * ones(1, NN); 

% subtract the neuropil
caCorrected = bsxfun(@times, neu, coefNeu);

if ops.recomputeKernel
    tlag                     = 1;
    [kernel, mtau]           = estimateKernel(ops, caCorrected, tlag);
    
%     [kernel, mtau, coefNeu]  = estimateKernel(ops, ca - coefNeu * neu, tlag);
    
    fprintf('Timescale determined is %4.4f samples \n', mtau);
end
kernel = normc(kernel(:));
%%
Params = [1 6 1 2e4]; % parameters of deconvolution
npad        = 250;    % padding

kernelS     = repmat(kernel, 1, NN);
dcell       = cell(NN,1);

Fsort       = sort(my_conv2(caCorrected, ceil(mtau), 1), 1, 'ascend');
baselines   = Fsort(ceil(NT/20), :); % estimate baseline

% determine and subtract the baselines
F1      = caCorrected - bsxfun(@times, ones(NT,1), baselines);

% normalize corrected traces
sd      = .5 * std(F1 - my_conv2(F1, 2, 1), [], 1); % ceil(ops.fs/2)
F1      = bsxfun(@rdivide, F1 , 1e-12 + sd);

% run the deconvolution to get fs etc
parfor icell = 1:size(ca,2)
    [~,dcell{icell}] = ...
        single_step_single_cell(F1(:,icell), Params, kernelS(:,icell), NT, npad,dcell{icell});
end
%     disp([mean(F1(:)>1e-12)  mean(err(:))   ])

sp = zeros(size(F1));
% rescale baseline contribution
for icell = 1:size(ca,2)
    dcell{icell}.c                      = dcell{icell}.c * sd(icell);
    dcell{icell}.neuropil_coefficient   = coefNeu(icell);
    dcell{icell}.baseline               = baselines(icell);
    
    sp(dcell{icell}.st) = dcell{icell}.c;
end

%%

