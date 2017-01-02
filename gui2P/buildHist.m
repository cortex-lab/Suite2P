function model = buildHist(prior, X, statLabels)

if ~isempty(prior)
    model = prior;
    
%     for j = 2:size(X,2)
%         grids(:,j-1) = linspace(myICDF(X(:,j), .02), myICDF(X(:,j), .98), 100);
%     end
    
    grids = prior.grids;
    if ~isempty(X)
        iscell  = X(:,1)>0.5;
        ix      = {iscell, ~iscell};
        
        for i = 1:numel(ix)
            ncells      = sum(ix{i});
            if ncells>0
                Xhist(:,:,i) = smoothDistro(X(ix{i}, 2:end), grids);
                
                f(i)            = ncells/(ncells + prior.n(i));
                Xhist(:,:,i)    = (1-f(i))*prior.hist(:,:,i) + f(i) *  Xhist(:,:,i);
            end
        end
        
        model.hist   = Xhist;
        model.grids  = grids;
    end
else
    if isempty(X)
        error('prior and database cannot both be empty')
    end
    
end
end

function Xhist = smoothDistro(X, grids)

Xhist = zeros(size(grids));
% smooth with sig proportional to pX
sig = 10;

for j = 1:size(grids,2)
    xbin = X(:,j);
    xbin(xbin<grids(1,j)) = grids(1,j);
    xbin(xbin>grids(end,j)) = grids(end,j);
    
    pX0 = histcounts(xbin, grids(:,j));
    pX = sig/2 * my_conv2(pX0, sig, 2);
    
    for k = 1:numel(pX0)
        L = zeros(size(grids,1), 1);
        if pX0(k)>0
            L(k, 1) = pX0(k);
            Xhist(:,j) = Xhist(:,j) + my_conv2(L, sig / pX(k), 1);
        end
    end
end
Xhist = bsxfun(@rdivide, Xhist, sum(Xhist,1));
end

function Xt = myICDF(X, p)

if p<0 || p>1
   error('quartile level must be between 0 and 1') 
end

[Xsort] = sort(X(:));

Xt = Xsort(max(1, ceil(p*numel(Xsort))));
end