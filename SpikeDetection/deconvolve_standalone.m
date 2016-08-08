function dcell = deconvolve_standalone(Ff, ops, Fneu)
% outputs a cell array of deconvolved spike times and amplitudes.

% you need to set: ops.toolbox_path, 

% the following options are set to their defaults
ops.imageRate    = getOr(ops, {'imageRate'}, 30); % total image rate (over all planes)
ops.sensorTau    = getOr(ops, {'sensorTau'}, 2); % approximate timescale in seconds
ops.sameKernel   = getOr(ops, {'sameKernel'}, 1); % 1 for same kernel per plane, 0 for individual kernels (not recommended)
ops.maxNeurop    = getOr(ops, {'maxNeurop'}, Inf);
ops.recomputeKernel    = getOr(ops, {'recomputeKernel'}, 0);


flag_neurop = 1;
if nargin<3
    flag_neurop = 0;
end

% load the initialization of the kernel    
load(fullfile(ops.toolbox_path, 'SpikeDetection\kernel.mat'));

% the basis functions should depend on timescale of sensor and imaging rate
mtau             = ops.imageRate * ops.sensorTau/ops.nplanes; 

coefNeu = .9; % initialize neuropil subtraction coef with 0.9
if flag_neurop
    sd0   = std(Ff-coefNeu*Fneu - my_conv2(Ff-coefNeu*Fneu, 2, 1), [], 1); % scale
    Fneu = 2 * Fneu ./ repmat(1e-5 + sd0, size(Ff,1), 1);
else
    sd0   = std(Ff - my_conv2(Ff, 2, 1), [], 1); % scale
end
Ff   = 2 * Ff ./ repmat(1e-5 + sd0, size(Ff,1), 1);

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
        [B(:,icell), err(icell)] = get_neurop(Ff(:,icell), F1(:, icell), kernelS(:,icell), npad, Fneu(:,icell));
    end
%     Ball{iter} = B;
%     save('C:\Users\Marius\Documents\MATLAB\Ball.mat', 'Ball')
    
    % determine and subtract the neuropil
    if flag_neurop
        F1 = Ff - bsxfun(@times, Fneu, B(3,:)) - bsxfun(@times, ones(NT,1), B(2,:));
        B(3, :) = min(maxNeurop, B(3,:));
    else
        F1 = Ff - bsxfun(@times, ones(NT,1), B(2,:));
    end
    
    % normalize signal by the sd of the neuropil-subtracted trace
    sd   = std(F1 - my_conv2(F1, 2, 1), [], 1);
    F1   = bsxfun(@rdivide, F1 , 1e-12 + sd/2);

    % run the deconvolution to get fs etc
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
    fprintf('fraction of non-zero bins %2.2f \n', mean(F1(:)>1e-12))    
    
    % multiply back by the sd
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
    
    % print the mean error
    fprintf('mean error %2.2f \n', mean(err(:))
   
end
%
fprintf('finished... \n')

% rescale baseline contribution
for icell = 1:size(Ff,2)
    dcell{icell}.B = B(:,icell);
    dcell{icell}.B(2) = dcell{icell}.B(2) * sd(icell);
end

%%

