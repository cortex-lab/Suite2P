function [ops, stat, res] = fastClustNeuropilCoef(ops, U, Sv)
% 
U =  reshape(U, [], size(U,ndims(U)));
iplane = ops.iplane;
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
Nkiter = round(linspace(ops.Nk0, ops.Nk, niter-2));
Nkiter(end+1:(niter+1)) = ops.Nk;

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
% addpath('C:\CODE\Github\Suite2P\SpikeDetection')
% f0      = 3; % sampling frequency
% load('D:\CODE\MariusBox\SpikeDetection\kernel.mat')

% Finv = (Fs*Fs')\Fs;

Uneu = U;
Ireg = diag([ones(Nk,1); zeros(nBasis,1)]);
%
tic
for k = 1:niter
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
    
    % recompute cell activities
    vs = zeros(nSVD, Nk); % initialize cell activities
    for i = 1:Nk
        ix = find(iclust==i);
        if numel(ix)>0
            vs(:,i) = M(ix) * Ucell(:, ix)';
        end
    end
    %
%     Ff = Fs' * vs;
%     vs = Finv * max(0, Ff);
%     [dcell, Ffr] = run_deconvolution2(Ff, f0, kernel);    
%     vs = Finv * Ffr;
    
    vs = bsxfun(@rdivide, vs, sum(vs.^2,1).^.5 + 1e-8);% normalize activity vectors
    
    % recompute pixels' assignments
    xs          = vs' * Ucell;
    [M, iclust] = max(xs,[],1);
    Uneu        = U - bsxfun(@times, M, vs(:,iclust)); %what's left over for neuropil model
    err(k)      = sum(sum((Uneu-neuropil).^2)).^.5;
    
    if 1
        %---------------------------------------------%
        xs(iclust + (0:Nk:numel(xs)-1)) = 0;
        [M2, iclust2] = max(abs(xs),[],1);
        
        dM = M - M2;
        indrem = Nkiter(k) - Nkiter(k+1);
        dMk = 0*ones(Nk,1);
        icl = cell(Nk,1);
        for j = 1:Nk
            icl{j} = iclust==j;
            if ~isempty(icl{j})
                dMk(j) = sum(dM(icl{j}));
            end
        end
        
        vlk = zeros(Nk, 1);
        vlk(~ison) = Inf;
        while indrem>0
            [Xmin, imin]          = min(dMk + vlk);
            if isinf(Xmin)
                Nkiter(k+1) = sum(ison);
                break;
            end
            newi               = iclust2(icl{imin});
            iclust(icl{imin})  = newi;
            M(icl{imin})       = M2(icl{imin});
            dMk(unique(newi))  = Inf;
            dMk(imin)          = Inf;
            dMk(unique(iclust(iclust2==imin))) ...
                = Inf;
            
            ison(imin) = 0;
            indrem             = indrem - 1;
        end
    end
    %---------------------------------------------%
    
    
    if (rem(k,10)==1 || k==niter) && ops.ShowCellMap
%         imagesc(reshape(PixL, Ly, Lx), [0 2])
%         drawnow
%         
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

%%
newindx = cumsum(ison);
iclust  = newindx(iclust);
Nk      = numel(unique(iclust));
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

if ~exist(ops.ResultsSavePath, 'dir')
    mkdir(ops.ResultsSavePath)
end
save(sprintf('%s/F_%s_%s_plane%d_Nk%d.mat', ops.ResultsSavePath, ...
    ops.mouse_name, ops.date, iplane, Nk),  'ops', 'res', 'stat')

%%0
% sk = skewness(F,[],2);
% [~, isk] = sort(sk, 'descend');
% clf
% for i = 1:20
%    plot(5*i + zscore(F(isk(i), :)))
%    hold all
% end
% axis tight
