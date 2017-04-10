function stat = run_deconvolution3(ops, dat)
% outputs a cell array of deconvolved spike times and amplitudes.

stat = dat.stat;

% construct Ff and Fneu
Ff = [];
Fneu = [];
for j = 1:numel(dat.Fcell)
    Ff   = cat(1, Ff, dat.Fcell{j}');
    Fneu = cat(1,Fneu, dat.FcellNeu{j}');
end


% the basis functions should depend on timescale of sensor and imaging rate
mtau             = ops.imageRate * ops.sensorTau/ops.nplanes; 

coefDefault = .9; % initialize neuropil subtraction coef with 0.8

Params = [1 1 .5 2e4]; %type of deconvolution, Th, Thi(nner loop), max Nspikes

% f0 = (mtau/2); % resample the initialization of the kernel to the right number of samples
kernel = exp(-[1:ceil(5*mtau)]'/mtau);
%
npad = 250;
[NT, NN] = size(Ff);

ops.fs            = getOr(ops, 'fs', ops.imageRate/ops.nplanes);

if ~isfield(stat, 'neuropilCoefficient')
    [coefNeu, inomax] = my_ica(Ff, Fneu, ops.fs, coefDefault, ops.maxNeurop);
else
   coefNeu = [stat.neuropilCoefficient]; 
end

Ff = Ff - bsxfun(@times, Fneu, coefNeu(:)');

% determine and subtract the running baseline
[F1, ~, sn] = estimate_baseline(Ff, ops);

if ops.recomputeKernel
    [kernel, mtau] = get1Dkernel(F1, ops.fs);

%     [kernel, mtau]  = estimateKernel(ops, Ff, tlag);
        
    fprintf('Timescale determined is %4.4f samples \n', mtau);
end
kernel = normc(kernel(:));
%%
kernelS   = repmat(kernel, 1, NN);
dcell = cell(NN,1);


switch ops.deconvType
    case 'OASIS'
        [sp, ca, coefs, B] = wrapperOASIS(ops, F, N);
    case 'L0'

end
% run the deconvolution to get fs etc
parfor icell = 1:size(Ff,2)
    [F1(:,icell),dcell{icell}] = ...
        single_step_single_cell(F1(:,icell), Params, kernelS(:,icell), NT, npad,dcell{icell});
end

fprintf('Percent bins with spikes %2.4f \n', 100*mean(F1(:)>0))



% rescale baseline contribution
for icell = 1:size(Ff,2)
    stat(icell).c                      = dcell{icell}.c;     % spike amplitudes
    stat(icell).st                     = dcell{icell}.st; % spike times    
    stat(icell).neuropilCoefficient    = coefNeu(icell); 
    stat(icell).noiseLevel             = sn(icell);     % noise level
    stat(icell).kernel                 = dcell{icell}.kernel;     % noise level
end


%%

