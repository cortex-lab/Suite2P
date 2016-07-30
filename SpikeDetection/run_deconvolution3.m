function [dcell, isroi] = run_deconvolution3(ops, dat)
% outputs a cell array of deconvolved spike times and amplitudes.
% Optionally output this in matrix form Ffr (very sparse).

% load the initialization of the kernel    
load(fullfile(ops.toolbox_path, 'SpikeDetection\kernel.mat'));

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
sd0   = std(Ff-coefNeu*Fneu - my_conv2(Ff-coefNeu*Fneu, 2, 1), [], 1);
Ff   = 2 * Ff ./ repmat(1e-5 + sd0, size(Ff,1), 1);
Fneu = 2 * Fneu ./ repmat(1e-5 + sd0, size(Ff,1), 1);

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
nt0 = numel(kernel);

taus = mtau * [1/8 1/4 1/2 1 2];
Nbasis = numel(taus);
kerns = zeros(nt0, Nbasis);
for i = 1:Nbasis
    kerns(:,i) = exp(-[0:nt0-1]/taus(i)) - exp(-[0:nt0-1]/(taus(i)/2));
end
kerns = normc(kerns);
%

kernelS = repmat(max(0, kernel), 1, NN);

maxNeurop = ops.maxNeurop;

err = zeros(NN, 1);

FfA = zeros(Nbasis, NN);
AtA = zeros(Nbasis, Nbasis,NN);
F1 = zeros(size(Ff), 'single');
B2 = zeros(Nbasis, NN);

dcell = cell(NN,1);

tic
for iter = 1:10
    plot(kernelS)
    drawnow
    
    fprintf('%2.2f sec, iter %d... ', toc, iter)
   
    % determine neuropil and cell contributions
    parfor icell = 1:size(Ff,2)
        [B(:,icell), err(icell)] = get_neurop(Ff(:,icell), Fneu(:,icell), F1(:, icell), kernelS(:,icell), npad);
    end
%     Ball{iter} = B;
%     save('C:\Users\Marius\Documents\MATLAB\Ball.mat', 'Ball')
    
    % determine and subtract the neuropil
    F1 = Ff - bsxfun(@times, Fneu, B(3,:)) - bsxfun(@times, ones(NT,1), B(2,:));
    
    % normalize signal 
    sd   = std(F1 - my_conv2(F1, 2, 1), [], 1);
    F1   = bsxfun(@rdivide, F1 , 1e-12 + sd/2);

    % run the deconvolution to get fs etc\
    if iter==10
        parfor icell = 1:size(Ff,2)
            [~, ~, ~,dcell{icell}] = ...
                single_step_single_cell(F1(:,icell), Params, kernelS(:,icell), kerns, NT, npad,dcell{icell});
        end
        break;
    else
        parfor icell = 1:size(Ff,2)
            [FfA(:,icell), AtA(:,:,icell), F1(:,icell)] = ...
                single_step_single_cell(F1(:,icell), Params, kernelS(:,icell), kerns, NT, npad);
        end
    end
    mean(F1(:)>1e-12)
    
    F1   = bsxfun(@times, F1 , 1e-12 + sd/2);
    
    % determine new kernel parameters
    if ops.recomputeKernel
        if ops.sameKernel
            B2 = mean(FfA,2)' /mean(AtA,3);
            B2 = repmat(B2', 1, NN);
        else
            for i = 1:NN
                B2(:,i) = FfA(:,i)'/ AtA(:,:,i);
            end
        end
        
        kernelS = normc(max(0, kerns * B2));
        kernelS = realign_kernels(kernelS, mtau/2 + ops.nplanes/ops.imageRate);
        kernelS = normc(max(0,kernelS));
    end
    
    mean(err(:))
   
end
%
fprintf('finished... \n')

% rescale baseline contribution
for icell = 1:size(Ff,2)
    dcell{icell}.B = B(:,icell);
    dcell{icell}.B(2) = dcell{icell}.B(2) * sd(icell);
end

%%

