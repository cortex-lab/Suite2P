function model = buildaLDA(prior, X, statLabels)

if ~isempty(prior)
    model = prior;
    
    if ~isempty(X)
        iFoot = strcmp('footprint', statLabels);
        mix = X(:,iFoot)>1;
        
        iscell  = X(:,1)>0.5;
        ix      = {iscell, ~iscell &  mix, ~iscell &  ~mix};
        
        mu = prior.mu;
        V  = prior.V;
        for i = 1:3
            ncells      = sum(ix{i});
            if ncells>0
                f(i)        = ncells/(ncells + prior.n(i));
                mu(i,:)     = (1-f(i))*mu(i,:)   + f(i) * mean(X(ix{i}, 2:end), 1);
                V(:,:,i)    = (1-f(i))*V(:,:,i)  + f(i) *  cov(X(ix{i}, 2:end));
            end
        end
        
        model.mu = mu;
        model.V  = V;
    end
else
    if isempty(X)
        error('prior and database cannot both be empty')
    end
    iFoot = strcmp('footprint', statLabels);
    mix = X(:,iFoot)>1;
    
    iscell  = X(:,1)>0.5;
    ix = {iscell, ~iscell &  mix, ~iscell &  ~mix};
    
    for i = 1:3
        mu(i,:)     = mean(X( ix{i}, 2:end), 1);
        V(:,:,i)    = cov(X( ix{i}, 2:end));
    end
    
    model.mu = mu;
    model.V  = V;
end
