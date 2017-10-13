function [sp, ca, sd, F1, tproc, kernel] = deconvolution_standalone(ops, ca, neu, gt)
% takes as input the calcium and (optionally) neuropil traces,  
% both NT by NN (number of neurons).

% specify in ops the following options, or leave empty for defaults
%       fs = sampling rate
%       sensorTau  = timescale of sensor, if recomputeKernel = 0
% additional options can be specified (mostly for linking with Suite2p, see below).

ops.imageRate    = getOr(ops, {'imageRate'}, 30); % total image rate (over all planes)
ops.nplanes      = getOr(ops, {'nplanes'}, 1); % how many planes at this total imaging rate
ops.sensorTau    = getOr(ops, {'sensorTau'}, 2); % approximate timescale in seconds
ops.maxNeurop    = getOr(ops, {'maxNeurop'}, Inf); % maximum allowed neuropil contamination coefficient. 

ops.lam              = getOr(ops, 'lam', 3);

if ops.lam<.25
    warning('L0 penalty with lambda < 1 is *very* slow to solve... \n')
elseif ops.lam<1
    warning('L0 penalty with lambda < 1 is slow to solve... \n')
end
if ops.lam<1e-10
    error('Penalty lambda is < 1e-10. Cannot solve. Typical values are 1 to 10... \n')
end

% the kernel should depend on timescale of sensor and imaging rate
ops.fs           = getOr(ops, 'fs', ops.imageRate/ops.nplanes);
mtau             = ops.fs * ops.sensorTau; 

if nargin<3 || isempty(neu)
    neu = zeros(size(ca));
end

Params = [1 ops.lam 1 2e4]; %parameters of deconvolution

% f0 = (mtau/2); % resample the initialization of the kernel to the right number of samples
kernel = exp(-[1:ceil(5*mtau)]'/mtau);
%
npad        = 250;
[NT, NN]    = size(ca);

%%

% determine and subtract the running baseline
if getOr(ops, 'runningBaseline', 0)
    Fbase    = ca;
    NT      = size(Fbase,1);
    ntBase  = 2*ceil(ops.running_baseline * ops.fs/2)+1;
    Fbase    = cat(1, Fbase((ntBase-1)/2:-1:1, :), Fbase, Fbase(end:-1:end-(ntBase-1)/2, :));
    Fbase   = my_conv2(Fbase, 3, 1); % medfilt1(Fneu(:), 5);
    Fbase   = ordfilt2(Fbase, 1, true(ntBase,1));
    
    Fbase = Fbase((ntBase-1)/2 + [1:NT], :);
    
%     y = Fneu - Fbase;    
elseif getOr(ops, 'doBaseline', 1)
    Fsort       = my_conv2(ca, ceil(ops.fs), 1);
    Fsort       = sort(Fsort, 1, 'ascend');
    baselines   = Fsort(ceil(NT/20), :);
    Fbase   = bsxfun(@times, ones(NT,1), baselines);
else
    Fbase = zeros(size(ca));
end
F1 = ca - Fbase;

% normalize by the noise. EXTRA FACTOR OF 1/2 compared with OASIS.
sd   = 1/2 * 1/sqrt(2) * std(F1(2:end, :) - F1(1:end-1, :), [], 1);
F1   = bsxfun(@rdivide, F1 , 1e-12 + sd);


%%

if ops.sensorTau>1e-3
    if isfield(ops, 'fracTau')
        kernel1 = exp(-[1:ceil(5*mtau)]'/mtau);
        kernel2 = exp(-[1:ceil(5*mtau)]'/(mtau * ops.fracTau));
        
        kernel = kernel1 - kernel2;
    end
    kernel = normc(kernel(:));
    kernel = repmat(kernel, 1, size(F1,2));
else
    kernel = zeros(ceil(5*ops.fs), size(F1,2));
    for k = 1:size(F1,2)
        g = estimate_time_constant(F1(:,k), 1);
        g = g(1);
        
        kernel(:, k) = g.^[1:ceil(5*ops.fs)];
    end
    kernel = normc(kernel);
end

if nargin>3
    % get GT kernel
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
    kernel = kerns * coefs;
    kernel = normc(kernel(:));
end


sp = zeros(size(F1), 'single');

tic
% run the deconvolution to get fs etc\
parfor icell = 1:size(ca,2)
    sp(:,icell) = ...
        single_step_single_cell(F1(:,icell), Params, kernel(:, icell), NT, npad);
end
tproc = toc;

parfor icell = 1:size(ca,2)
    ca(:,icell) = filter(kernel(:,icell), 1, sp(:,icell));    
end
ca   = bsxfun(@times, ca, sd);
F1   = bsxfun(@times, F1, sd);



%%

