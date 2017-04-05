function [ops, stat, model] = sourcery2(ops, U, model)

ops.fig         = getOr(ops, 'fig', 1);
ops.ThScaling   = getOr(ops, 'ThScaling', 1);

U0 =  reshape(U, [], size(U,ndims(U)))';

Ly = numel(ops.yrange);
Lx = numel(ops.xrange);

[nSVD, Npix] = size(U0);

U0 = reshape(U0, nSVD, Ly, Lx);

S = getNeuropilBasis(ops, Ly, Lx, 'raisedcosyne'); % 'raisedcosyne', 'Fourier'
S = normc(S);
nBasis = size(S,2);

%
d0   = ceil(ops.diameter); % expected cell diameter in pixels

sig = ceil(d0/4); 

% data Ucell is nMaps by Ly by Lx
dx = repmat([-d0:d0], 2*d0+1, 1);
dy = dx';

rs = dx.^2 + dy.^2 - d0^2;
dx = dx(rs<=0);
dy = dy(rs<=0);

mPix    = zeros(numel(dx), 1e4);
mLam    = zeros(numel(dx), 1e4);
mLam0   = zeros(numel(dx), 1e4);

iter = 0;
icell = 0;
r = [];

L   = sparse(Ly*Lx, 0);
LtU = zeros(0, nSVD);
LtS = zeros(0, nBasis);

% neu     = StS\StU;
% Ucell   =  U0 - reshape(neu' * S', size(U0));

U0 = U0(:,:);
Ucell = U0;

nBasis = size(S,2);
PixL = ones(1, Lx * Ly)';
Sm = bsxfun(@times, S, PixL);

Th = 2 * ops.ThScaling;
%% 
while 1
    iter = iter + 1;    
    %re-estimate neuropil
    StS = Sm' * Sm;
    StU = Sm' * Ucell';
    Lam = (StS + 1e-3 * eye(nBasis)) \ StU;
    
    % recompute neuropil pixel contribution
    neuropil = Lam' * S';
    if 1
        PixL = mean(bsxfun(@times, neuropil, Ucell), 1);
        PixL = bsxfun(@rdivide, PixL, mean(neuropil.^2,1));%         PixL = max(0, PixL);
    else
        PixL  = PixL';
    end
    neuropil = bsxfun(@times, neuropil, PixL);
    
    Ucell = Ucell - gather(neuropil);
     
    PixL  = PixL';
    Sm = bsxfun(@times, S, PixL);
    
    % re-estimate masks
    L   = sparse(Ly*Lx, icell);
    for j = 1:icell
        ipos = find(mPix(:,j)>0);
        ipix = mPix(ipos,j);
        
        Usub = Ucell(:, ipix)+ codes(j, :)' * mLam(ipos,j)';
        
        lam = max(0, codes(j, :) * Usub);
%         threshold pixels
        lam(lam<max(lam)/5) = 0;
        
        mLam(ipos,j) = lam;
        
%         extract biggest connected region of lam only
        mLam(:,j) = normc(getConnected(mLam(:,j), rs));
        lam = mLam(ipos,j);
        
%         normalize and multiply by movie SD for display purposes
        mLam0(ipos,j) = lam .* model.sdmov(ipix);
        
        L(ipix,j) = lam;
        
        LtU(j, :) = U0(:,ipix) * lam;
        LtS(j, :) = lam' * Sm(ipix,:);
        
        Ucell(:, ipix) = Usub - (Usub * lam)* lam';
    end
  
    % what variance is left in total?
    err(iter) = 1e3*mean(Ucell(:).^2);
    fprintf('%d total ROIs, err %4.4f, thresh %4.4f \n', icell, err(iter), Th)
    
    % residual is smoothed at every iteration
    us = my_conv2_circ(reshape(Ucell, nSVD, Ly, Lx), sig, [2 3]);
    
    % compute log variance at each location
    V = double(sq(mean(us.^2,1)));    

    if 1
        if iter==1
            um = sq(mean(Ucell.^2,1));
            um = my_conv2_circ(reshape(um, Ly, Lx), sig, [1 2]);
        end
        V = V./double(um);
    end
    
    % do the morpholrogical opening trick
    if iter==1
        lbound = -my_min2(-my_min2(V, d0), d0);
        lbound = my_conv2(lbound, d0, [1 2]);
    end
    
    V = V - lbound;
    
    if iter==1
       % find indices of all maxima  in plus minus 1 range
        maxV    = -my_min(-V, 1, [1 2]);
        ix      = (V > maxV-1e-10);
        
        % threshold is the mean peak, times a potential scaling factor
        pks = V(ix);

        Th  = ops.ThScaling * median(pks(pks>1e-4)); 
    end
    
    % find local maxima in a +- d neighborhood
    maxV = -my_min(-V, 2*d0, [1 2]);
    
    % find indices of these maxima above a threshold
    ix  = (V > maxV-1e-10) & (V > Th);
    ind = find(ix);
    
    if iter==1
       Nfirst = numel(ind); 
       V0 = V;
       ops.Vcorr = V0;       
    end
    
    % just in case this goes above original value
    V = min(V, V0);
    
    if numel(ind)==0 || numel(ind)<Nfirst * getOr(ops, 'stopSourcery', 1/20)
        break;
    end
    
    new_codes = normc(us(:, ind));
    
    ncells = icell;
    LtU(ncells+size(new_codes,2), nSVD) = 0;
    
    % each source needs to be iteratively subtracted off
    for i = 1:size(new_codes,2)
        icell = icell + 1;
        [ipix, ipos] = getIpix(ind(i), dx, dy, Lx, Ly);
        
        Usub = Ucell(:, ipix);
        
        lam = max(0, new_codes(:, i)' * Usub);        
        
        % threshold pixels
        lam(lam<max(lam)/5) = 0;
         
        mPix(ipos,icell) = ipix;
        mLam(ipos,icell) = lam;
        
        % extract biggest connected region of lam only
        mLam(:,icell)   = normc(getConnected(mLam(:,icell), rs));
        lam             = mLam(ipos,icell) ;
        
        L(ipix,icell) = lam;
        
        LtU(icell, :) = U0(:,ipix) * lam;
        LtS(icell, :) = lam' * Sm(ipix,:);
        
        Ucell(:, ipix) = Usub - (Usub * lam)* lam';
    end    
    
   % ADD NEUROPIL INTO BIG REGRESSION HERE        
    StS     = Sm' * Sm;
    StU     = Sm' * U0';
    LtL     = full(L'*L);
    
    codes   = ([LtL LtS; LtS' StS]+ 1e-3 * eye(icell+nBasis)) \ [LtU; StU];
    neu     = codes(icell+1:end,:);    
    codes   = codes(1:icell,:);    
    
    Vnew    = sq(sum(Ucell.^2,1));
    Ucell   = U0 - reshape(double(codes') * L', nSVD, Ly*Lx);    
    
end
%%
if ops.fig
    figure
    subplot(1,2, 1);
    imagesc(V0, [0 4*Th])
    axis off
    
    subplot(1,2, 2);
    imagesc(V, [0 Th])
    axis off    
    
    drawnow
end

fprintf('%d total ROIs, err %4.4f, thresh %4.4f \n', icell, err(end), Th)

mLam  =  mLam(:, 1:icell);
mLam0 = mLam0(:, 1:icell);
mPix  =  mPix(:, 1:icell);
%%

% subtract off neuropil only
Ucell = U0 - reshape(neu' * Sm', size(U0));

% populate stat with cell locations and footprint
stat = getFootprint(ops, codes, Ucell, mPix, mLam, mLam0);

% compute compactness of ROIs
stat = anatomize(ops, mPix, mLam, stat);

figure
[~, iclust, lam] = drawClusters(ops, r, mPix, mLam, Ly, Lx);

model.L     = L;
model.S     = S;
model.LtS   = LtS;
model.LtL   = LtL;
model.StS   = StS;

% get anatomical projection weights
stat = weightsMeanImage(ops, stat, model);

