function dcell = deconvolution_standalone(ops, ca, neu)
% takes as input the calcium and (optionally) neuropil traces,  
% both NT by NN (number of neurons).

% outputs a cell array dcell containing spike times (dcell.st) and amplitudes
% (dcell.c). dcell.B(3) is the neuropil contamination coefficient. 
% dcell.B(2) is an estimate of the baseline. 

ops.imageRate    = getOr(ops0, {'imageRate'}, 30); % total image rate (over all planes)
ops.nplanes      = getOr(ops0, {'nplanes'}, 1); % how many planes at this total imaging rate
ops.sensorTau    = getOr(ops0, {'sensorTau'}, 2); % approximate timescale in seconds
ops.sameKernel   = getOr(ops0, {'sameKernel'}, 1); % 1 for same kernel per plane, 0 for individual kernels (does not work right now)
ops.maxNeurop    = getOr(ops0, {'maxNeurop'}, Inf); % maximum allowed neuropil contamination coefficient. 
ops.recomputeKernel = getOr(ops0, {'recomputeKernel'}, 1); % whether to estimate kernel from data

% the kernel should depend on timescale of sensor and imaging rate
ops.fs           = getOr(ops, 'fs', ops.imageRate/ops.nplanes);
mtau             = ops.fs * ops.sensorTau; 

if nargin<3 || isempty(neu)
    neu = zeros(size(ca));
end

Params = [1 3 1 2e4]; %type of deconvolution, Th, Thi(nner loop), max Nspikes

% f0 = (mtau/2); % resample the initialization of the kernel to the right number of samples
kernel = exp([-1:ceil(5*mtau)]'/mtau);
%
npad        = 250;
[NT, NN]    = size(ca);
coefNeu = .8 * ones(1,NN); % initialize neuropil subtraction coef with 0.8

caCorrected = ca - bsxfun(@times, neu, coefNeu);

if ops.recomputeKernel
    tlag                     = 1;
    [kernel, mtau]           = estimateKernel(ops, caCorrected, tlag);
    
%     [kernel, mtau, coefNeu]  = estimateKernel(ops, ca - coefNeu * neu, tlag);
    
    fprintf('Timescale determined is %4.4f samples \n', mtau);
end
kernel = normc(kernel(:));
%%
B = zeros(3, NN);

kernelS   = repmat(kernel, 1, NN);
dcell = cell(NN,1);

tic


% determine and subtract the neuropil
F1 = caCorrected - bsxfun(@times, ones(NT,1), baselines);

% normalize signal
sd   = std(F1 - my_conv2(F1, 2, 1), [], 1);
F1   = bsxfun(@rdivide, F1 , 1e-12 + sd/2);

% run the deconvolution to get fs etc\
parfor icell = 1:size(ca,2)
    [sp(:,icell),dcell{icell}] = ...
        single_step_single_cell(F1(:,icell), Params, kernelS(:,icell), NT, npad,dcell{icell});
end

%     disp([mean(F1(:)>1e-12)  mean(err(:))   ])

F1   = bsxfun(@times, F1 , 1e-12 + sd/2);

% rescale baseline contribution
for icell = 1:size(ca,2)
    dcell{icell}.c    = dcell{icell}.c * sd(icell)/2;
    dcell{icell}.B    = B(:,icell);
    dcell{icell}.B(2) = dcell{icell}.B(2) * sd(icell);
end

%%

