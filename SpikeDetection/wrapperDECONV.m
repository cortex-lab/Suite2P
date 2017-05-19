function [sp, ca, coefs, B, sd, ops] = wrapperDECONV(ops, F, N)

sd = std(F,1,1);

N = bsxfun(@rdivide, N, sd);
F = bsxfun(@rdivide, F, sd);

F(:, isnan(sum(N,1))) = nanmean(F,2) * ones(1,sum(isnan(sum(F,1))));
N(:, isnan(sum(N,1))) = nanmean(N,2) * ones(1,sum(isnan(sum(N,1))));

NN      = size(F,2);
B       = zeros(2, NN);

coefs   = .7 * ones(1, NN);
if ops.estimateNeuropil
    niter = 5;
else
    niter = 1;
end

% maximum neuropil coefficient
ops.maxNeurop = getOr(ops, 'maxNeurop', 1);

for k = 1:niter
    Fsub = bsxfun(@minus, F, bsxfun(@times, N, coefs));
    
    % at each iteration run the deconvolution on the residual
    switch getOr(ops, 'deconvType', 'L0')
        case 'OASIS'
            [sp, ca, sd2] = OASISpreprocess(ops,  Fsub);
        case 'L0'
            [sp, ca, sd2] = deconvolution_standalone(ops, Fsub);
        otherwise
            error('Please choose deconvolution method as OASIS or L0')
    end
    
    [NT, NN] = size(F);
    clear B
    
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

sd = sd.*sd2;

ca = bsxfun(@times, ca, sd);
