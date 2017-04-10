function [sp, ca, sd] = OASISpreprocess(ops, ca, neu)
% takes as input the calcium and (optionally) neuropil traces,  
% both NT by NN (number of neurons).

% specify in ops the following options, or leave empty for defaults
%       fs = sampling rate
%       sensorTau  = timescale of sensor, if recomputeKernel = 0
% additional options can be specified (mostly for linking with Suite2p, see below).

ops.imageRate    = getOr(ops, {'imageRate'}, 30); % total image rate (over all planes)
ops.nplanes      = getOr(ops, {'nplanes'}, 1); % how many planes at this total imaging rate
ops.sensorTau    = getOr(ops, {'sensorTau'}, 2); % approximate timescale in seconds
ops.sameKernel   = getOr(ops, {'sameKernel'}, 1); % 1 for same kernel per plane, 0 for individual kernels (does not work right now)
ops.maxNeurop    = getOr(ops, {'maxNeurop'}, Inf); % maximum allowed neuropil contamination coefficient. 
ops.recomputeKernel = getOr(ops, {'recomputeKernel'}, 1); % whether to estimate kernel from data

% the kernel should depend on timescale of sensor and imaging rate
ops.fs           = getOr(ops, 'fs', ops.imageRate/ops.nplanes);

if nargin<3 || isempty(neu)
    neu = zeros(size(ca));
end

[NT NN] = size(ca);

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
else
    Fsort       = my_conv2(ca, ceil(ops.fs), 1);
    Fsort       = sort(Fsort, 1, 'ascend');
    baselines   = Fsort(ceil(NT/20), :);
    Fbase   = bsxfun(@times, ones(NT,1), baselines);
end

F1 = ca - Fbase;

% normalize signal
sd   = 1/sqrt(2) * std(F1(2:end, :) - F1(1:end-1, :), [], 1);
F1   = bsxfun(@rdivide, F1 , 1e-12 + sd);

sp = zeros(size(F1), 'single');
ca = zeros(size(F1), 'single');

lam         = getOr(ops, 'lam', 3);
mtau        = ops.fs * ops.sensorTau; 
g           = exp(-1/mtau);
parfor n = 1:size(F1,2)
    [ca(:,n), sp(:,n)] = oasisAR1(F1(:,n)', g, lam);
%     [c_oasis, sp(:,n)] = oasisAR2(F1(:,n)', [], lam);
%     [c_oasis, sp(:,n)] = constrained_oasisAR1(F1(:,n)', g);
    %         [c_oasis, sp(:,n)] = oasisAR1(Raw{id}(:,n)', g);
end

ca   = bsxfun(@times, ca, sd);

%%

