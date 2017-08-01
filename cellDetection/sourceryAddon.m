% this should be run after sourcery to refine the spatial extent of the
% masks (removes the effect of smoothing on the raw data)

% reshape U to be (nMaps x Y x X)
U0 =  reshape(U2, [], size(U2,ndims(U2)))';
Ly = numel(ops.yrange);
Lx = numel(ops.xrange);
[nSVD, Npix] = size(U0);
U0 = reshape(U0, nSVD, Ly, Lx);

% regress maps onto basis functions and subtract neuropil contribution
StU     = S'*U0(:,:)'; % covariance of neuropil with spatial masks
neu     = StS\StU;

nBasis = size(S,2);

% set to 0 the masks, to be re-estimated
mLam    = zeros(numel(dx), 1e4);
L       = sparse(Ly*Lx, icell); 

for iter = 1:3
    % subtract off everything
    Ucell = U0 - reshape(neu' * S', size(U0)) - reshape(double(codes') * L', size(U0));    
    
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
        
        Ucell(:, ipix) = Usub - (Usub * lam)* lam';
    end
    err(iter) = mean(Ucell(:).^2);
    
    % ADD NEUROPIL INTO REGRESSION HERE    
    Ucell = Ucell + reshape(neu' * S', size(U0));    
    
    StU     = S'*Ucell(:,:)'; % covariance of neuropil with spatial masks
    neu     = StS\StU;
    
    fprintf('%d total ROIs, err %4.4f, thresh %4.4f \n', icell, err(iter), Th)
    if ops.fig           
        figure(3)
        [~, iclust, lam] = drawClusters(ops, r, mPix, mLam, Ly, Lx);
        
        drawnow
    end
end
%%