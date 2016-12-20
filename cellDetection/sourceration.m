function [] = sourcery(U, Sv, Ly, Lx, ops)

% clearvars -except U Sv Ly Lx ops
U0 =  reshape(U, [], size(U,ndims(U)));
U0 = bsxfun(@times, U0, Sv'.^.5)';

Ly = numel(ops.yrange);
Lx = numel(ops.xrange);

[nSVD, Npix] = size(U0);

U0 = reshape(U0, nSVD, Ly, Lx);

% U0 = U0 - my_conv2(U0, 5, [2 3]);

ops.TileFactor = 4;
S = getNeuropilBasis(ops, Ly, Lx, 'raisedcosyne'); % 'raisedcosyne', 'Fourier'
S = normc(S);
nBasis = size(S,2);

StU = S'*U0(:,:)';
StS = S'*S;
%%
d0   = 12; % expected cell diameter

sig = d0/4; % 
d   = 2*d0; % max expected cell diameter
% Th  = 3; %1e-6

% data Ucell is nMaps by Ly by Lx

dx = repmat([-d/2:d/2], d+1, 1);
dy = dx';

rs = dx.^2 + dy.^2;
mask = exp(-rs/(2*sig^2));
mask = mask/mean(mask(:));
dx = dx(rs<=(d/2)^2);
dy = dy(rs<=(d/2)^2);


mPix = zeros(numel(dx), 1e4);
mLam = zeros(numel(dx), 1e4);


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
    us = my_conv2(Ucell, sig, [2 3]);
    
    % compute log variance at each location
    V = sq(mean(us.^2,1));
    V = double(V);
    
    um = sq(mean(Ucell.^2,1));
    um = my_conv2(um, sig, [1 2]);
    
    V = V./um ;
%     V = log(V./um);
    
     V = double(V);
    % do the morpholrogical opening trick
    if iter==1
        lbound = -my_min(-my_min(V, d/2, [1 2]), d/2, [1 2]);
        V = V - lbound;
        
        % find indices of all maxima  in plus minus 1 range
        maxV = -my_min(-V, 1, [1 2]);
        ix  = (V > maxV-1e-10); 
            
        % threshold is the mean peak, times a potential scaling factor
        pks = V(ix);
        Th = mean(pks);        
    else
        V = V - lbound;
    end
    
    % find local maxima in a +- d neighborhood
    maxV = -my_min(-V, d, [1 2]);
    
    % find indices of these maxima above a threshold
    ix  = (V > maxV-1e-10) & (V > Th);
    ind = find(ix);
    
    if numel(ind)==0
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
        
        lam = normc(lam(:));
         
        mLam(ipos,j) = lam;
        
        L(ipix,j) = lam;
        
        LtU(j, :) = U0(:,ipix) * lam;
        LtS(j, :) = lam' * S(ipix,:);
        
        Ucell(:, ipix) = Usub - (Usub * lam)* lam';
    end
    
    err(iter) = mean(Ucell(:).^2*1e6);
    fprintf('%d total cells %d new, err %4.4f \n', icell, size(new_codes,2), err(iter))
    
    Vnew = sq(sum(Ucell.^2,1));
    figure(1)
    subplot(1,2,1)
    imagesc(V, [Th 2*Th])
    
    subplot(1,2,2)
    imagesc(Vnew)
    
    figure(2)
    r = drawClusters(r, mPix, mLam, Ly, Lx);
    
    drawnow
end
%%
xlx         = repmat(-ceil(d0):1:ceil(d0), 2*ceil(d0)+1, 1);
rgrid       = sqrt(xlx.^2 + xlx'.^2);
rgridsort   = sort(rgrid(:), 'ascend');

ds = (dx.^2 + dy.^2).^.5;

for j = 1:icell
    ipos = find(mPix(:,j)>0);
    ipix = mPix(ipos,j);
    lam  = mLam(ipos,j);
        
    gpix = lam>max(lam)/4;
    rd(j) = mean(ds(gpix));
    rd0(j) = mean(rgridsort(1:sum(gpix)));    
end

ThL = d0/4;
ThH = d0/2;
iscell = rd<ThH & rd>ThL;
ThL = 1.0;
ThH = 1.1;
iscell = iscell & rd./rd0<ThH & rd./rd0>ThL;

figure(3)
drawClusters(r, mPix(:, iscell), mLam(:, iscell), Ly, Lx);
title('final')
%%





