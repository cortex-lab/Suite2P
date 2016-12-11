function [dcell, isroi] = run_deconvolution3(ops, dat)
% outputs a cell array of deconvolved spike times and amplitudes.
% Optionally output this in matrix form Ffr (very sparse).

% load the initialization of the kernel    
% load(fullfile(ops.toolbox_path, 'SpikeDetection\kernel.mat'));

if isfield(dat.stat, 'igood')
   isroi = logical([dat.stat.igood]); 
else
    if isfield(dat, 'cl') && isfield(dat.cl, 'iscell')
        isroi = dat.cl.iscell;
    else
        isroi = [dat.stat.mrs]./[dat.stat.mrs0]<dat.clustrules.Compact & ...
            [dat.stat.npix]>dat.clustrules.MinNpix & [dat.stat.npix]<dat.clustrules.MaxNpix;
    end
    isroi = logical(isroi);
end

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

coefDefault = .9; % initialize neuropil subtraction coef with 0.8

Params = [1 1 .5 2e4]; %type of deconvolution, Th, Thi(nner loop), max Nspikes

% f0 = (mtau/2); % resample the initialization of the kernel to the right number of samples
kernel = exp([-1:ceil(5*mtau)]'/mtau);
%
npad = 250;
[NT, NN] = size(Ff);

ops.fs            = getOr(ops, 'fs', ops.imageRate/ops.nplanes);
[coefNeu, inomax] = my_ica(Ff, Fneu, ops.fs, coefDefault);

Ff = Ff - bsxfun(@times, Fneu, coefNeu(:)');

% determine and subtract the running baseline
[F1, ~] = estimate_baseline(Ff, ops);

if ops.recomputeKernel
    [kernel, mtau] = get1Dkernel(F1, ops.fs);

%     [kernel, mtau]  = estimateKernel(ops, Ff, tlag);
        
    fprintf('Timescale determined is %4.4f samples \n', mtau);
end
kernel = normc(kernel(:));
%%
kernelS   = repmat(kernel, 1, NN);
dcell = cell(NN,1);


% run the deconvolution to get fs etc\
parfor icell = 1:size(Ff,2)
    [F1(:,icell),dcell{icell}] = ...
        single_step_single_cell(F1(:,icell), Params, kernelS(:,icell), NT, npad,dcell{icell});
end

mean(F1(:)>0)

% rescale baseline contribution
for icell = 1:size(Ff,2)
    dcell{icell}.c                      = dcell{icell}.c;    
    dcell{icell}.neuropil_coefficient   = coefNeu(icell);
end


%%

