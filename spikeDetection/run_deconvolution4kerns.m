function [fs, ca, B, Bk, kernel] = run_deconvolution4(ops, Ff, Fneu)
% this version (4) assumes all ROIs should be deconvolved

% outputs an array filled with deconvolved spike times and amplitudes, 
% as well as smoothed traces


if nargin<3 || isempty(Fneu)
    Fneu= zeros(size(Ff));
end


% initialize the kernel with 1 second timescales
ops.sensorTau = getOr(ops, 'sensorTau', 2);
ops.imageRate = getOr(ops, 'imageRate', 30);
ops.nplanes   = getOr(ops, 'nplanes', 1);
ops.maxNeurop = getOr(ops, 'maxNeurop', 1);

ops.fs        = getOr(ops, 'fs', ops.imageRate / ops.nplanes);

% timescale in frames
mtau    =  ops.fs * ops.sensorTau;

if isfield(ops, 'running_baseline') && ops.running_baseline>0
    Fneu = Ff;
    NT = size(Ff,1);
    ntBase = 2*ceil(ops.running_baseline * ops.fs/2)+1;
    Fneu = cat(1, Fneu((ntBase-1)/2:-1:1, :), Fneu, Fneu(end:-1:end-(ntBase-1)/2, :));
    Fbase = my_conv2(Fneu, 3, 1); % medfilt1(Fneu(:), 5);
    Fbase = ordfilt2(Fbase, 1, true(ntBase,1));

    Fneu = Fneu((ntBase-1)/2 + [1:NT], :);
    Fbase = Fbase((ntBase-1)/2 + [1:NT], :);
    
    Ff = Fneu - Fbase;
end


% kernel
kernel  = exp(-[1:ceil(5*mtau)]/mtau);
kernel  = kernel(:);

% kernel = cat(1, zeros(ceil(5*mtau), 1), kernel, zeros(ceil(5*mtau), 1));
% kernel = kernel - expfilt(kernel', ceil(2*mtau))';
% kernel = kernel((ceil(5*mtau) + 1):end);

kernel  = normc(kernel(:));

kernelS = repmat(kernel, 1, size(Ff,2));

coefNeu = .8; % initialize neuropil subtraction coef with 0.8
Params  = [1 6 1 2e4]; %type of deconvolution, Th, Thi(nner loop), max Nspikes


nt0 = numel(kernel);

taus = mtau * [1/4 1];
Nbasis = numel(taus);
kerns = zeros(nt0, Nbasis);
for i = 1:Nbasis
    kerns(:,i) = exp(-[0:nt0-1]/taus(i)) - exp(-[0:nt0-1]/(taus(i)/2));
end
kerns = normc(kerns);


% Ff = Ff  - expfilt(Ff', ceil(2*mtau))';

npad        = 250; %padding with zeros
[NT, NN]    = size(Ff);

%%
B   = zeros(3, NN);

err = zeros(NN, 1);

ca = zeros(size(Ff));
fs = zeros(size(Ff));

Fsort = sort(my_conv2(Fneu, ceil(mtau), 1), 1, 'ascend');
B(2,:) = Fsort(ceil(NT/20), :);


tic
for iter = 1:10
%     fprintf('%2.2f sec, iter %d... ', toc, iter)
   
    % subtract the neuropil
    F1 = Ff - bsxfun(@times, Fneu, B(3,:)) - bsxfun(@times, ones(NT,1), B(2,:));
    
    % normalize signal 
    sd   = .5 * std(F1 - my_conv2(F1, 2, 1), [], 1);
    F1   = bsxfun(@rdivide, F1 , 1e-12 + sd);

    % run the deconvolution to get fs etc
    parfor icell = 1:size(Ff,2)
        [fs(:,icell), ca(:,icell)] = ...
            single_step_single_cell2(F1(:,icell), Params, kernelS(:,icell), npad);
    end
    
    Fk = zeros(NT+nt0-1, size(fs,2), Nbasis);
    for k = 1:Nbasis
        kk = kerns(:,k);
        parfor j = 1:size(fs,2)    
            Fk(:,j, k) = conv(kk, fs(:,j));
        end
    end
    Fk = Fk(1:NT, :, :);
    Fk = reshape(Fk, [], Nbasis);
    Bk = (Fk'*Fk)\(Fk' * F1(:));
    kernel = normc(kerns * Bk);
    kernelS = repmat(kernel, 1, size(Ff,2));
    
    
    plot(kernel)
    drawnow
    
    ca   = bsxfun(@times, ca , 1e-12 + sd);
    fs   = bsxfun(@times, fs , 1e-12 + sd);
    
   % determine neuropil and cell contributions
    parfor icell = 1:size(Ff,2)
        [~, err(icell)] = get_neurop(Ff(:,icell)-ca(:, icell), Fneu(:,icell));
    end
     
    err = err + Params(2) * ((sd').^2).* mean(fs>1e-10, 1)';
    
%     B(2,:) = max(0, B(2,:));
    [mean(fs(:)>1e-12) mean(err(:)) mean(B(2,:))]
end
%
% fprintf('finished... \n')

ca = bsxfun(@plus, ca, B(2,:));

%%

