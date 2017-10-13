function [sp, ca, sd, F1, Fbase] = OASISpreprocess(ops, ca, neu, gt)
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

[NT , ~] = size(ca);

% determine and subtract the running baseline
Fbase = zeros(size(ca), 'single');

if getOr(ops, 'runningBaseline', 0)   
    Fbase    = ca;
    if getOr(ops, 'zBaseline', 0)    
        [~, ix] = sort(ops.zdrift);
        Fbase = Fbase(ix, :);
    end
    
    ntBase  = 2*ceil(ops.runningBaseline * ops.fs/2)+1;
    Fbase   = cat(1, Fbase((ntBase-1)/2:-1:1, :), Fbase, Fbase(end:-1:end-(ntBase-1)/2, :));
        
    Fbase   = my_conv2(Fbase, 1, 1); 
    Fbase   = movmin(Fbase, ntBase,1);
    Fbase   = movmax(Fbase, ntBase,1);
    
    Fbase = Fbase((ntBase-1)/2 + [1:NT], :);    
    
    if getOr(ops, 'zBaseline', 0)
        Fbase(ix,:) = Fbase;     
    end
end
F1          = ca - Fbase;

% Fsort       = my_conv2(F1, ceil(ops.fs), 1);
% Fsort       = sort(Fsort, 1, 'ascend');
% baselines   = Fsort(ceil(NT/20), :);
% Fbase2      = bsxfun(@times, ones(NT,1), baselines);
% F1          = F1 - Fbase2;
% Fbase       = Fbase + Fbase2; 

% F1         = F1./max(std(F1,1,1), Fbase);

% normalize signal
sd   = 1/sqrt(2) * std(F1(2:end, :) - F1(1:end-1, :), [], 1);
F1   = bsxfun(@rdivide, F1 , 1e-12 + sd);

sp = zeros(size(F1), 'single');
ca = zeros(size(F1), 'single');

lam         = getOr(ops, 'lam', 0);
mtau        = ops.fs * ops.sensorTau; 
g           = exp(-1/mtau);

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
    
    kerns = exp(-bsxfun(@rdivide, [1:ceil(5*mtau)]',  mtau * exp(linspace(-3, 1.5, 11))));
    err = Inf * ones(size(kerns,2));
       
    for j = 1:size(kerns,2)-1
        for i = j+1:size(kerns,2)
            Y = kerns(:, [i j]);
            B = (Y'*Y)\(Y'*kernel);
            err(j,i) = mean((kernel - Y*B).^2);
        end
    end
    [i, j] = find(err==min(err(:)));
    
    Y = kerns(:, [i j]);
    B = (Y'*Y)\(Y'*kernel);
    plot(Y*B)
    hold all
    plot(kernel)
    hold off
    drawnow
    
    g = kerns(1, [j i]);
        
    g = [(g(1)+g(2))  -g(1)*g(2)];
    parfor n = 1:size(F1,2)
        [ca(:,n), sp(:,n)] = oasisAR2(F1(:,n)', g, lam);
    end

elseif isfield(ops, 'fracTau') 
    g           = [exp(-1/mtau) exp(-1/(ops.fracTau*mtau))];
    g = [(g(1)+g(2))  -g(1)*g(2)];
    parfor n = 1:size(F1,2)
        [ca(:,n), sp(:,n)] = oasisAR2(F1(:,n)', g, lam);
    end
elseif g>1e-3
    parfor n = 1:size(F1,2)
        [ca(:,n), sp(:,n)] = oasisAR1(F1(:,n)', g, lam);
        %     [c_oasis, sp(:,n)] = constrained_oasisAR1(F1(:,n)', g);
        %         [c_oasis, sp(:,n)] = oasisAR1(Raw{id}(:,n)', g);
    end
else 
    for n = 1:size(F1,2)
        [ca(:,n), sp(:,n)] = oasisAR1(F1(:,n)', [], lam);
        %     [c_oasis, sp(:,n)] = constrained_oasisAR1(F1(:,n)', g);
        %         [c_oasis, sp(:,n)] = oasisAR1(Raw{id}(:,n)', g);
    end  
end

ca   = bsxfun(@times, ca, sd);
sp   = bsxfun(@times, sp, sd);
F1   = bsxfun(@times, F1, sd);

if getOr(ops, 'zBaseline', 0)
%     ca   = ca.*Fbase;
    ca = ca + Fbase;
end

end

% spfilt0 = ones([size(kerns,2) size(F1)]);
% for j = 1:size(kerns,2)
%     spfilt0(j, :,:) = filter(kerns(:,j), 1, gt);
% end
% 
% err = Inf * ones(size(kerns,2));
% for j = 1:size(kerns,2)-1
%     for i = j+1:size(kerns,2)
%         spfilt = spfilt0([j i], :, :);
%         
%         sts = zeros(size(spfilt,1));
%         stF = zeros(size(spfilt,1), 1);
%         for k = 1:size(spfilt,3)
%             sts = sts + spfilt(:,:,k) * spfilt(:,:,k)';
%             stF = stF + spfilt(:,:,k) * F1(:,k);
%         end
%         coefs   = sts \ stF;
%         pred    = coefs' * spfilt(:,:);
%         err(j,i) = mean((pred(:) - F1(:)).^2);
%     end
% end
% 
% [i, j] = find(err==min(err(:)));