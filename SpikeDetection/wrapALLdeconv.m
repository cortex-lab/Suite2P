function [ops, stat, dat] = wrapALLdeconv(ops, stat, Fcell, FcellNeu)

% set up options for deconvolution
ops.imageRate    = getOr(ops0, {'imageRate'}, 30); % total image rate (over all planes)
ops.sensorTau    = getOr(ops0, {'sensorTau'}, 2); % approximate timescale in seconds
ops.sameKernel   = getOr(ops0, {'sameKernel'}, 1); % 1 for same kernel per plane, 0 for individual kernels (not recommended)
ops.sameKernel   = getOr(ops0, {'sameKernel'}, 1);
ops.maxNeurop    = getOr(ops0, {'maxNeurop'}, Inf); % maximum neuropil coefficient (default no max)
ops.recomputeKernel    = getOr(ops0, {'recomputeKernel'}, 0);
ops.deconvNeuropil = getOr(ops0, {'deconvNeuropil'}, 0);   % whether to deconvolve the neuropil as well

fprintf('Spike deconvolution, plane %d... \n', iplane)

% construct Ff and Fneu
Ff = [];
Fneu = [];
for j = 1:numel(Fcell)
    Ff   = cat(1, Ff, Fcell{j}');
    Fneu = cat(1,Fneu, FcellNeu{j}');
end

ops.fs                  = ops.imageRate/ops.nplanes;
ops.estimateNeuropil    = getOr(ops0, 'estimateNeuropil', 1);
ops.runningBaseline     = 0;

[sp, ~, coefs,~, sd, ops] = wrapperDECONV(ops, Ff, Fneu);

if ops.deconvNeuropil
    ops.estimateNeuropil = 0;
    spNeu = wrapperDECONV(ops, Fneu);
end

for j = 1:size(Ff,2)
    stat(j).neuropilCoefficient = coefs(j);
    stat(j).noiseLevel          = sd(j);
end

nCum = 0;
for j = 1:length(dat.Fcell)
    dat.sp{j} = sp(nCum + [1:size(dat.Fcell{j},2)], :)';
    if ops.deconvNeuropil
        dat.spNeu{j} = spNeu(nCum + [1:size(dat.Fcell{j},2)], :)';
    end
    nCum = nCum + size(dat.Fcell{j},2);
end