function  [Ypred, ps] = evaluateHist(model, X)

grids = model.grids;
ntypes = size(model.hist,3);

logps = zeros(size(X,1), ntypes);
for j = 1:size(model.grids,2)
    xp = X(:, j+1);
    xp(xp<grids(1,j))   = grids(1,j);
    xp(xp>grids(end,j)) = grids(end,j);
    
    [~, ~, ibin]    = histcounts(xp, grids(:,j));
    
    logps           = logps + log(sq(model.hist(ibin, j, :)));
end

ps = ones(1,ntypes)/ntypes; 
for j = 1:10
    L    = bsxfun(@plus, logps, log(ps));
    L    = bsxfun(@minus, L, max(L, [], 2));
    rs   = exp(L) + 1e-5;
    rs   = bsxfun(@rdivide, rs, sum(rs,2));
    ps   = mean(rs,1);
end

Ypred = rs(:,1);    
