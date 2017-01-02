function [Ypred, ps] = evaluateLDA(model, X)

ntypes = size(model.mu,1);
logps = zeros(size(X,1), ntypes);

for i = 1:ntypes
    Xm          = bsxfun(@minus, X(:, 2:end), model.mu(i,:));
    logps(:,i)  = logps(:,i)- sum((Xm/model.V(:,:,i)) .* Xm, 2)/2 - ...
        sum(log(abs(eigs(model.V(:,:,i)))))/2;
end

%% optimize the mixing probabilities
ps = ones(1,ntypes)/ntypes; %model.p;
for i = 1:10
   L    = bsxfun(@plus, logps, log(ps));
   L    = bsxfun(@minus, L, max(L, [], 2));
   rs   = exp(L) + 1e-5;
   rs   = bsxfun(@rdivide, rs, sum(rs,2));
   ps   = mean(rs,1);
end
% ps
Ypred = rs(:,1);    

