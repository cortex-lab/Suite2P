function [sp, ca, coefs, B, sd, ops, baselines] = wrapperDECONV(ops, F, N)

% F contains the NT by NN somatic traces
% optional: N contains the NT by NN neuropil traces

% IMPORTANT: specify in ops the following options
%       ops.fs = sampling rate
%       ops.sensorTau  = timescale of sensor, if recomputeKernel = 0
%       ops.estimateNeuropil = 0 or 1 (default is 1, if N is given as input)

% outputs:
%     sp: the deconvolved traces
%     ca: the reconstructed calcium
%     coefs: the neuropil coefficients, if inferred
%     B,sd,ops: for Suite2p use
%     baselines: inferred baselines

sd = std(F,1,1);
F = bsxfun(@rdivide, F, sd);
[NT, NN] = size(F);

% maximum neuropil coefficient
ops.maxNeurop = getOr(ops, 'maxNeurop', 1);

runBaseSave = getOr(ops, 'runningBaseline', 0);

if nargin<=2
   ops.estimateNeuropil = 0;
end
if getOr(ops, 'estimateNeuropil', 0)
    niter = 5;    
    ops.runningBaseline = 0;    
else
    niter = 1;
end

F(:, isnan(sum(F,1))) = nanmean(F,2) * ones(1,sum(isnan(sum(F,1))));    
if nargin>2
    N = bsxfun(@rdivide, N, sd);
    N(:, isnan(sum(F,1))) = nanmean(N,2) * ones(1,sum(isnan(sum(F,1))));    
else
    niter = 1;
end

coefs   = min(.8, ops.maxNeurop) * ones(1, NN);
% coefs(:) = .7;
B       = zeros(2, NN);


for k = 1:niter
    if nargin>2
        Fsub = bsxfun(@minus, F, bsxfun(@times, N, coefs));
    else
        Fsub = F;
    end
    
    % at each iteration run the deconvolution on the residual
    switch getOr(ops, 'deconvType', 'OASIS')
        case 'OASIS'
            if exist('oasisAR1.m', 'file') ~= 2
               warning('Could not find oasisAR1.m in your path. Please download the separate OASIS github repository and add to your path.') ;
               warning('Refer to instructions at top of example master file for more help.')
               error('OASIS not found.')
            end
            [sp, ca, sd2, Fbase] = OASISpreprocess(ops,  Fsub);
        case 'L0'
            [sp, ca, sd2] = deconvolution_standalone(ops, Fsub);
        otherwise
            error('Please choose deconvolution method as OASIS or L0')
    end
    
    if niter>1 && k<niter
        Y = F - ca;
        
        % given the spike trace, infer the baseline and the neuropil contribution
        for j = 1:NN
            % baseline should be order one, just like N and F
            X = cat(2, ones(NT,1), N(:,j));
            B(:, j) = ((X'*X)/size(X,1) + 1e-3 * eye(size(X,2))) \ ...
                ((X'*Y(:,j))/size(X,1));
        end
        
        B(isnan(B)) = 1;
        coefs       = min(B(2,:), ops.maxNeurop);
        coefs       = max(0, coefs);
       
        disp([median(coefs) std(coefs(:))])
    end
end

% rescale the calcium and deconvolved
ca = bsxfun(@times, ca, sd);
sp = bsxfun(@times, sp, sd);
baselines = bsxfun(@times, mean(Fbase,1), sd);

sd = sd.*sd2;

ops.runningBaseline = runBaseSave;
