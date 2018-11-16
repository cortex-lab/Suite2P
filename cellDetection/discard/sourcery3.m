% function [ops, stat, model] = sourcery2(ops, U, model)

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

%%
while 1
    iter = iter + 1;    
   
    %recompute neuropil    
    Sm = bsxfun(@times, S, PixL);
    StU = Sm' * Ucell';
    StS = Sm' * Sm;
    Lam = (StS + 1e-4 * eye(nBasis)) \ StU;
    
    % recompute neuropil pixel contribution
    neuropil = Lam' * S';
    PixL = mean(bsxfun(@times, neuropil, Ucell), 1);
    PixL = bsxfun(@rdivide, PixL, mean(neuropil.^2,1));
%     PixL = max(0, PixL);
    %     PixL  = PixL';
    neuropil = bsxfun(@times, neuropil, PixL);
    PixL  = PixL';
    
            Ucell = Ucell - gather(neuropil);
            Uneu  = U0    - gather(neuropil);
%     Uneu = Ucell - gather(neuropil);
%     err(iter) = mean(Uneu(:).^2);
%     disp(err(iter))
%     
%     figure(2)
%     imagesc(reshape(mean(Uneu.^2,1), Ly, Lx))
%     
%     figure(3)
%     imagesc(reshape(PixL, Ly, Lx))
%     
%     drawnow
    
    % re-estimate masks
    L   = sparse(Ly*Lx, icell);
    for j = 1:icell
        ipos = find(mPix(:,j)>0);
        ipix = mPix(ipos,j);
        
        Usub = Ucell(:, ipix)+ codes(j, :)' * mLam(ipos,j)';
        
        lam = max(0, codes(j, :) * Usub);
        % threshold pixels
        lam(lam<max(lam)/5) = 0;
        
        mLam(ipos,j) = lam;

        % extract biggest connected region of lam only
        mLam(:,j) = normc(getConnected(mLam(:,j), rs));
        lam = mLam(ipos,j);
        
        L(ipix,j) = lam;
        
        LtU(j, :) = Uneu(:,ipix) * lam;
        
        Ucell(:, ipix) = Usub - (Usub * lam)* lam';
        
        % normalize and multiply by movie SD for display purposes
        mLam0(ipos,j) = lam .* model.sdmov(ipix);
    end
    
    % residual is smoothed at every iteration
    us = my_conv2_circ(reshape(Ucell, nSVD, Ly, Lx), sig, [2 3]);
    
    % compute log variance at each location
    V = sq(mean(us.^2,1));
    V = double(V);
    
    um = sq(mean(Ucell.^2,1));
    um = my_conv2_circ(reshape(um, Ly, Lx), sig, [1 2]);
    
    V = V./um ;
%     V = log(V./um);
    
     V = double(V);
    % do the morpholrogical opening trick
    if iter==1
        lbound = -my_min2(-my_min2(V, d0), d0);
    end
    
    V = V - lbound;
    
    if iter==1        
        % find indices of all maxima  in plus minus 1 range
        maxV    = -my_min(-V, 1, [1 2]);
        ix      = (V > maxV-1e-10);
        
        % threshold is the mean peak, times a potential scaling factor
        pks = V(ix);

        Th  = ops.ThScaling * median(pks(pks>1e-4));
        
        V0 = V;
        ops.Vcorr = V0;
    end
    
    % find local maxima in a +- d neighborhood
    maxV = -my_min(-V, 2*d0, [1 2]);
    
    % find indices of these maxima above a threshold
    ix  = (V > maxV-1e-10) & (V > Th);
    ind = find(ix);
    
    if iter==1
       Nfirst = numel(ind); 
    end
    
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
        
        LtU(icell, :) = Uneu(:,ipix) * lam;
        LtS(icell, :) = lam' * S(ipix,:);
    end   
    %
   % ADD NEUROPIL INTO BIG REGRESSION HERE    
    LtL     = full(L'*L);
    
    codes = (LtL + 1e-3 * eye(icell))\LtU;
    
    % subtract off everything
    Ucell = Uneu - reshape(double(codes') * L', size(U0));    
   
    %
    err(iter) = mean(Ucell(:).^2);
    
    Vnew = sq(sum(Ucell.^2,1));
    
    Ucell = U0 - reshape(double(codes') * L', nSVD, Ly*Lx);    
    
    fprintf('%d total ROIs, err %4.4f, thresh %4.4f \n', icell, err(iter), Th)
    
    if ops.fig
        figure(1)
        subplot(1,2, 1);
        imagesc(V0, [0 2*Th])
        axis off
        
        subplot(1,2, 2);
        imagesc(V, [0 2*Th])
        axis off
        
        figure(2)
        [~, iclust, lam] = drawClusters(ops, r, mPix, mLam, Ly, Lx);

        
        figure(3)
        imagesc(reshape(PixL, Ly, Lx))
        
        drawnow
    end
    
end
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
Ucell = U0 - reshape(neu' * S', size(U0));

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

