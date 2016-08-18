function [ops, stat, res] = fastClustNeuropilCoef2(ops, U, Sv, Vsvd,flag_save)
% 
if nargin<5
    flag_save = 0;  % don't save by default
end

U =  reshape(U, [], size(U,ndims(U)));
U = bsxfun(@times, U, Sv'.^.5)';
[nSVD, Npix] = size(U);
% Fs = bsxfun(@rdivide, Fs, Sv.^.5);

%% clustering options
Ly = numel(ops.yrange);
Lx = numel(ops.xrange);


ops.Nk0 = floor(sqrt(ops.Nk0))^2;

Nk      = ops.Nk0;
nsqrt = round(sqrt(Nk));

xs = repmat(round(linspace(1, nsqrt, Lx)), Ly, 1);
ys = repmat(round(linspace(1, nsqrt, Ly))', 1, Lx);
iclust = xs + (ys-1) * nsqrt;
iclust = iclust(:);

clear xs ys

niter   = ops.niterclustering;

% xs = repmat(1:Lx, Ly, 1);
% ys = repmat((1:Ly)', 1, Lx);
% 
% randx = rand(1, Nk) * Lx;
% randy = rand(1, Nk) * Ly;
% 
% dx = repmat(xs(:), 1, Nk) - repmat(randx, numel(xs(:)), 1);
% dy = repmat(ys(:), 1, Nk) - repmat(randy, numel(ys(:)), 1);
% 
% dxy = dx.^2 + dy.^2;
% [~, iclust] = min(dxy, [], 2);
%%
% clear dx dy

if ops.ShowCellMap
    figure( 'Units', 'pixels', 'position', [100 100 900 900])
    colormap('hsv')
    axes('position', [.05 .05 .925 .925])
    set(gcf, 'Color', 'w')
end

r   = rand(1, Nk);
Sat = ones(Ly, Lx);

err = zeros(niter,1);
ops.meanV = gather(sum(Sv)/(Ly*Lx));

Nk = ops.Nk0;
Nkiter = round(linspace(ops.Nk0, ops.Nk, niter+1));
% Nkiter(end+1:(niter+1)) = ops.Nk;

%%
Npix = Ly * Lx;
M = .0001 * ones(1, Npix, 'single');

ison = true(Nk,1);
TileFactor = getOr(ops, {'TileFactor'}, 1); % this option can be overwritten by the user
nTiles = ceil(TileFactor * (Ly+Lx)/2 / (10 * ops.diameter)); % neuropil is modelled as nTiles by nTiles 

xc = linspace(1, Lx, nTiles);
yc = linspace(1, Ly, nTiles);
yc = yc';
xs = 1:Lx;
ys = 1:Ly;

sigx = 4*(Lx - 1)/nTiles;
sigy = 4*(Ly - 1)/nTiles;

S = zeros(Ly, Lx, nTiles, nTiles, 'single');
for kx = 1:nTiles
    for ky = 1:nTiles
        cosx = 1+cos(2*pi*(xs - xc(kx))/sigx);
        cosy = 1+cos(2*pi*(ys - yc(ky))/sigy);
        cosx(abs(xs-xc(kx))>sigx/2) = 0;
        cosy(abs(ys-yc(ky))>sigy/2) = 0;
        
        S(:, :,ky, kx) = cosy' * cosx;
    end
end
S = reshape(S, [], nTiles^2);
S = normc(S);

nBasis = nTiles^2 ;
PixL = ones(1, Lx * Ly)';
%%
Uneu = U;
Ireg = diag([ones(Nk,1); zeros(nBasis,1)]);
%
flag_exit = 0;

nFeat = 10;

tic
k = 0;
while k<niter
    k = k+1;
    % recompute neuropil
    Sm = bsxfun(@times, S, PixL);
    StS = Sm' * Sm;
    StU = Sm' * Uneu';
    Lam = (StS + 1e-4 * eye(nBasis)) \ StU;
    
    % recompute neuropil pixel contribution 
    neuropil = Lam' * S';
    PixL = mean(bsxfun(@times, neuropil, Uneu), 1);
    PixL = bsxfun(@rdivide, PixL, mean(neuropil.^2,1));
    PixL = max(0, PixL);
    neuropil = bsxfun(@times, neuropil, PixL);
    Ucell = U - neuropil; %what's left over for cell model
    PixL = PixL';
    
    icl         = get_cl(iclust, Nk);
    iSortFeat   = get_isortfeat(icl, Nk, Ly, nFeat);
    
    % recompute cell activities
    vs = zeros(nSVD, Nk); % initialize cell activities
    for i = 1:Nk
        ix = icl{i}; %find(iclust==i);
        if numel(ix)>0
            vs(:,i) = M(ix) * Ucell(:, ix)';
        end
    end
    
    vs = bsxfun(@rdivide, vs, sum(vs.^2,1).^.5 + 1e-8);% normalize activity vectors
    
    
%     xs = zeros(nFeat, Npix);
%     xall = zeros(nFeat, Nk);
%     for i = 1:Nk
%         xs(:, icl{i}) = vs(:, iSortFeat(:,i))' * Ucell(:, icl{i});
%         xall(:,i) = sum(xs(:, icl{i}), 2);
%     end
%     
%     [M, ifeat] = max(xs,[],1);
%     iclust = iSortFeat(ifeat + nFeat * (iclust'-1));
    
    % recompute pixels' assignments
    xs          = vs' * Ucell;
    [M, iclust] = max(xs,[],1);
    Uneu        = U - bsxfun(@times, M, vs(:,iclust)); %what's left over for neuropil model
    err(k)      = sum(sum((Uneu-neuropil).^2)).^.5;
    
    % estimate merging matrix
    xall = zeros(nFeat, Nk);
    for i = 1:Nk        
        xall(:,i) = sum(xs(iSortFeat(:,i), icl{i}), 2);
    end
    if 1
        % find neighboring clusters
        icl         = get_cl(iclust, Nk);
        iSortFeat   = get_isortfeat(icl, Nk, Ly, nFeat);
        iSortFeat   = iSortFeat(2:nFeat, :);
        
        ccNeigh = bsxfun(@minus, xall(2:end, :), xall(1, :));
%         ccNeigh = ccNeigh./max(0, xall(2:end, :));

%         dF = Vsvd * vs;
%         CC = corrcoef(dF);        
%         ccNeigh = zeros(nFeat-1, Nk);
%         for i = 1:Nk
%             ccNeigh(:,i) = CC(iSortFeat(:, i), i);
%         end
 
        
        figure(2)
        hold off
        plot(sort(max(ccNeigh, [], 1)))
        
        hold on
        for j = 1:Nkiter(k)-Nkiter(k+1)
            [cMax, imax] = max(ccNeigh, [], 1);
            [cMaxMax, jmax] = max(cMax);
            if cMaxMax < 1/2 * median(cMax)
                if j==1
                    flag_exit = 1;
                end
                break;
            end
            imax = imax(jmax);
            imax = iSortFeat(imax, jmax);
            
            
            %             keyboard;
            iclust(icl{imax}) = jmax;
            ccNeigh(:, [imax jmax])    = -Inf;
            ccNeigh(iSortFeat(:)==imax) = -Inf;
            ccNeigh(iSortFeat(:)==jmax) = -Inf;
            
            Nk = Nk - 1;
            
        end
        plot(sort(max(ccNeigh, [], 1)))
        drawnow
        %         ibad = find(iclust>Nkiter(k+1));
        %         iclust(ibad) = ceil(rand(1,numel(ibad))*Nkiter(k+1));
        %         % find most correlated neighbor
        %
        %         % merge top N pairs
        %
        [~, ~, iclust] = unique(iclust);
%         iclust = iclust(:)';
        if k==niter && flag_exit==0
            k = k-1;
        end
    end
    %---------------------------------------------%
    
    
    if (rem(k,10)==1 || k==niter) && ops.ShowCellMap
%         imagesc(reshape(PixL, Ly, Lx), [0 2])
%         drawnow
%         
        figure(1)
        lam = M;
        for i = 1:Nk
            ix = find(iclust==i);
            nT0 = numel(ix);
            if nT0>0
                vM = lam(ix);
%                 vM = vM/sum(vM.^2)^.5;
                lam(ix) = vM;
            end
        end
%         V = max(0, min(10 * reshape(lam, Ly, Lx), 1));
        V = max(0, min(.5 * reshape(lam, Ly, Lx)/mean(lam(:)), 1));
        H = reshape(r(iclust), Ly, Lx);
        rgb_image = hsv2rgb(cat(3, H, Sat, V));
        imagesc(rgb_image)
        axis off
        drawnow
        fprintf('residual variance is %2.6f time %2.2f \n', err(k), toc)
    end
   
end

lam = M;
for i = 1:Nk
    ix = find(iclust==i);
    
    nT0 = numel(ix);
    if nT0>0
        vM = lam(ix);
        lam(ix) = vM/sum(vM.^2)^.5;
    end
end

%
clear res

res.iclust  = iclust;
res.M       = M;
res.S       = S;
res.lambda  = lam;

%
res.Ly  = Ly;
res.Lx  = Lx;
stat    = get_stat(res);


if flag_save
    if ~exist(ops.ResultsSavePath, 'dir')
        mkdir(ops.ResultsSavePath)
    end
    save(sprintf('%s/F_%s_%s_plane%d_Nk%d.mat', ops.ResultsSavePath, ...
        ops.mouse_name, ops.date, ops.iplane, Nk),  'ops', 'res', 'stat')
end

%%0
% sk = skewness(F,[],2);
% [~, isk] = sort(sk, 'descend');
% clf
% for i = 1:20
%    plot(5*i + zscore(F(isk(i), :)))
%    hold all
% end
% axis tight
