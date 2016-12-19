function [mPix, mLam, mLam0, model, Foot] = sourcery(ops, U, Sv)

% clearvars -except U Sv Ly Lx ops
U0 =  reshape(U, [], size(U,ndims(U)));
U0 = bsxfun(@times, U0, Sv'.^.5)';

Ly = numel(ops.yrange);
Lx = numel(ops.xrange);

[nSVD, Npix] = size(U0);

U0 = reshape(U0, nSVD, Ly, Lx);

% U0 = U0 - my_conv2(U0, 5, [2 3]);

ops.TileFactor = 4;
S = getNeuropilBasis(ops, Ly, Lx, 'Fourier'); % 'raisedcosyne', 'Fourier'
S = normc(S);
nBasis = size(S,2);

StU = S'*U0(:,:)';
StS = S'*S;
%%
d0   = ceil(ops.diameter); % expected cell diameter in pixels

sig = ceil(d0/4); 

% data Ucell is nMaps by Ly by Lx
dx = repmat([-d0:d0], 2*d0+1, 1);
dy = dx';

rs = dx.^2 + dy.^2;
dx = dx(rs<=d0^2);
dy = dy(rs<=d0^2);


mPix    = zeros(numel(dx), 1e4);
mLam    = zeros(numel(dx), 1e4);
mLam0   = zeros(numel(dx), 1e4);


iter = 0;
icell = 0;
r = [];

L   = sparse(Ly*Lx, 0);
LtU = zeros(0, nSVD);
LtS = zeros(0, nBasis);

neu     = StS\StU;
Ucell   =  U0 - reshape(neu' * S', size(U0));

%%
while 1
    iter = iter + 1;    
   
    % residual is smoothed at every iteration
    us = my_conv2_circ(Ucell, sig, [2 3]);
    
    % compute log variance at each location
    V = sq(mean(us.^2,1));
    V = double(V);
    
    um = sq(mean(Ucell.^2,1));
    um = my_conv2_circ(um, sig, [1 2]);
    
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
        Th  = median(pks(pks>1e-4));
        
        if ops.fig
            figure(1)
            my_subplot(1,3, 1, [.9 .9]);
            imagesc(V, [0 2*Th])
        end
    end
    
    % find local maxima in a +- d neighborhood
    maxV = -my_min(-V, 2*d0, [1 2]);
    
    % find indices of these maxima above a threshold
    ix  = (V > maxV-1e-10) & (V > Th);
    ind = find(ix);
    
    if iter==1
       Nfirst = numel(ind); 
    end
    
    if numel(ind)==0 || numel(ind)<Nfirst/20
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
        
        lam = normc(lam(:));
         
        mPix(ipos,icell) = ipix;
        mLam(ipos,icell) = lam;
        
        L(ipix,icell) = lam;
        
        LtU(icell, :) = U0(:,ipix) * lam;
        LtS(icell, :) = lam' * S(ipix,:);
    end    
    %
    % ADD NEUROPIL INTO REGRESSION HERE    
    LtL     = full(L'*L);
    codes = ([LtL LtS; LtS' StS]+ 0e-3 * eye(icell+nBasis))\[LtU; StU];
    neu   = codes(icell+1:end,:);    
    codes = codes(1:icell,:);
%     codes = (LtL+ 1e-3 * eye(icell))\LtU;    
    
    % subtract off everything
    Ucell = U0 - reshape(neu' * S', size(U0)) - reshape(double(codes') * L', size(U0));    
    
    L   = sparse(Ly*Lx, icell);
    for j = 1:icell        
        ipos = find(mPix(:,j)>0);
        ipix = mPix(ipos,j);        
        
        Usub = Ucell(:, ipix)+ codes(j, :)' * mLam(ipos,j)';
        
        lam = max(0, codes(j, :) * Usub);
        % threshold pixels
        lam(lam<max(lam)/5) = 0;
        
        
        mLam0(ipos,j) = lam;
        
        lam = normc(lam(:));
        mLam(ipos,j) = lam;
        
        L(ipix,j) = lam;
        
        LtU(j, :) = U0(:,ipix) * lam;
        LtS(j, :) = lam' * S(ipix,:);
        
        Ucell(:, ipix) = Usub - (Usub * lam)* lam';
    end
    
    err(iter) = mean(Ucell(:).^2*1e6);
    fprintf('%d total cells %d new, err %4.4f, thresh %4.4f \n', icell, size(new_codes,2), err(iter), Th)
    
    Vnew = sq(sum(Ucell.^2,1));
    
    if ops.fig
        figure(1)
        my_subplot(1,3, 3, [.9 .9]);
        imagesc(V, [Th 2*Th])
        axis off
        
        my_subplot(1,3, 2, [.9 .9]);
        imagesc(V, [0 2*Th])
        axis off
        
        figure(2)
        r = drawClusters(r, mPix, mLam, Ly, Lx);
        
        drawnow
    end
end

mLam = mLam(:, 1:icell);
mLam0 = mLam0(:, 1:icell);
mPix = mPix(:, 1:icell);
%%

% subtract off neuropil only
Ucell = U0 - reshape(neu' * S', size(U0));

%%
dx = repmat([-2*d0:2*d0], 4*d0+1, 1);
dy = dx';

rs = (dx.^2 + dy.^2).^.5;
dx = dx(rs<=2*d0);
dy = dy(rs<=2*d0);
rs = rs(rs<=2*d0);

frac = [0.15 0.25 0.33 0.5 .75];
Foot = zeros(icell, numel(frac));

% find maximum contamination distance for each ROI
for j = 1:icell
    ipos = find(mPix(:,j)>0);
    ipix = mPix(ipos,j);
    y0 = ceil(median(rem(ipix-1, Ly) + 1));
    x0 = ceil(median(ceil(ipix/Ly)));
    
    ivalid = find((x0 + dx)>=1 & (x0 + dx)<=Lx & (y0 + dy)>=1 & (y0 + dy)<=Ly);
    
    ipix = (y0+dy(ivalid)) + (x0 + dx(ivalid)-1) * Ly;
    proj = codes(j,:) * Ucell(:, ipix);
    
    for l = 1:numel(frac)
        Foot(j,l) = mean(rs(ivalid(proj>max(proj) * frac(l))));        
    end
end
Foot = Foot/ops.diameter;

model.L = L;
model.S = S;
model.LtS = LtS;
model.LtL = LtL;
model.StS = StS;



