function [sp, ca, coefs, B, sd, ops] = wrapperDECONV(ops, F, N)

if getOr(ops, 'estimateNeuropil', 0)
    niter = 5;
    % maximum neuropil coefficient
    ops.maxNeurop = getOr(ops, 'maxNeurop', 1);
else
    niter = 1;
end

sd = std(F,1,1);
F = bsxfun(@rdivide, F, sd);
[NT, NN] = size(F);

if nargin>2
    N = bsxfun(@rdivide, N, sd);
    F(:, isnan(sum(N,1))) = nanmean(F,2) * ones(1,sum(isnan(sum(F,1))));
    N(:, isnan(sum(N,1))) = nanmean(N,2) * ones(1,sum(isnan(sum(N,1))));
else
    niter = 1;
end

coefs   = .7 * ones(1, NN);
B       = zeros(2, NN);


for k = 1:niter
    if nargin>2
        Fsub = bsxfun(@minus, F, bsxfun(@times, N, coefs));
    else
        Fsub = F;
    end
    
    % at each iteration run the deconvolution on the residual
    switch getOr(ops, 'deconvType', 'L0')
        case 'OASIS'
            [sp, ca, sd2] = OASISpreprocess(ops,  Fsub);
        case 'L0'
            [sp, ca, sd2] = deconvolution_standalone(ops, Fsub);
        otherwise
            error('Please choose deconvolution method as OASIS or L0')
    end
    
    if niter>1
        Y = F - ca;
        
        % given the spike trace, infer the baseline and the neuropil contribution
        for j = 1:NN
            X = cat(2, ones(NT,1)/1e2, N(:,j));
            B(:, j) = ((X'*X)/size(X,1) + 1e-3 * eye(size(X,2))) \ ...
                ((X'*Y(:,j))/size(X,1));
        end
        
        B(isnan(B)) = 1;
        coefs       = min(B(2,:), ops.maxNeurop);
    end
end

% rescale the calcium and deconvolved
sd = sd.*sd2;
ca = bsxfun(@times, ca, sd);
