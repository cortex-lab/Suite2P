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

coefNeu = .9; % initialize neuropil subtraction coef with 0.8
sd0     = std(Ff-coefNeu*Fneu - my_conv2(Ff-coefNeu*Fneu, 2, 1), [], 1);
Ff      = 2 * Ff ./ repmat(1e-5 + sd0, size(Ff,1), 1);
Fneu    = 2 * Fneu ./ repmat(1e-5 + sd0, size(Ff,1), 1);

Params = [1 3 1 2e4]; %type of deconvolution, Th, Thi(nner loop), max Nspikes

% f0 = (mtau/2); % resample the initialization of the kernel to the right number of samples
kernel = exp([-1:ceil(5*mtau)]'/mtau);
%
npad = 250;
[NT, NN] = size(Ff);

if ops.recomputeKernel
    ops.fs = getOr(ops, 'fs', ops.imageRate/ops.nplanes);
    tlag                     = 1;
    [kernel, mtau]           = estimateKernel(ops, Ff - coefNeu * Fneu, tlag);
    
%     [kernel, mtau, coefNeu]  = estimateKernel(ops, Ff - coefNeu * Fneu, tlag);
    
    fprintf('Timescale determined is %4.4f samples \n', mtau);
end
kernel = normc(kernel(:));
%%
B = zeros(3, NN);

kernelS   = repmat(kernel, 1, NN);

err = zeros(NN, 1);

F1 = zeros(size(Ff), 'single');

dcell = cell(NN,1);

tic
niter = 10;
for iter = 1:niter
%     plot(kernelS)
%     drawnow
    
%     fprintf('%2.2f sec, iter %d... ', toc, iter)
   
    % determine neuropil and cell contributions
    parfor icell = 1:size(Ff,2)
        [B(:,icell), err(icell)] = ...
            get_neurop(Ff(:,icell), Fneu(:,icell), F1(:, icell), kernelS(:,icell), npad);
    end
   
    % determine and subtract the neuropil
    F1 = Ff - bsxfun(@times, Fneu, B(3,:)) - bsxfun(@times, ones(NT,1), B(2,:));
    
    % normalize signal 
    sd   = std(F1 - my_conv2(F1, 2, 1), [], 1);
    F1   = bsxfun(@rdivide, F1 , 1e-12 + sd/2);

    % run the deconvolution to get fs etc\
    if iter==niter
        parfor icell = 1:size(Ff,2)
            [F1(:,icell),dcell{icell}] = ...
                single_step_single_cell(F1(:,icell), Params, kernelS(:,icell), NT, npad,dcell{icell});
        end
    else
        parfor icell = 1:size(Ff,2)
            [F1(:,icell)] = ...
                single_step_single_cell(F1(:,icell), Params, kernelS(:,icell), NT, npad);
        end
    end
%     disp([mean(F1(:)>1e-12)  mean(err(:))   ])
    
    F1   = bsxfun(@times, F1 , 1e-12 + sd/2);   
end
%
% fprintf('finished... \n')

% rescale baseline contribution
for icell = 1:size(Ff,2)
    dcell{icell}.c    = dcell{icell}.c * sd(icell)/2;
    dcell{icell}.B    = B(:,icell);
    dcell{icell}.B(2) = dcell{icell}.B(2) * sd(icell);
end

%%

