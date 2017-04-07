% function [ops, stat, res, PixL] = fastClustNeuropilCoef(ops, U, Sv)

U =  reshape(U, [], size(U,ndims(U)));

U = bsxfun(@times, U, Sv'.^.5)';
[nSVD, Npix] = size(U);

%% clustering options
Ly = numel(ops.yrange);
Lx = numel(ops.xrange);


ops.Nk0 = floor(sqrt(ops.Nk0))^2;

Nk      = ops.Nk0;
niter   = ops.niterclustering;

%%
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

Nkiter = round(linspace(ops.Nk0, ops.Nk, niter-2));
Nkiter(end+1:(niter+1)) = ops.Nk;

%%
Npix = Ly * Lx;
M = ones(1, Npix, 'single');

ison = true(Nk,1);
S = getNeuropilBasis(ops, Ly, Lx, 'raisedcosyne'); % 'raisedcosyne', 'Fourier'

nBasis = size(S,2);
PixL = ones(1, Lx * Ly)';
%%

Uneu = U;
% Ireg = diag([ones(Nk,1); zeros(nBasis,1)]);


tic
for k = 1:niter
    %recompute neuropil    
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

%     [A Sv B] = svdecon(gpuArray(single(Uneu)));
%     neuropil = A(:, 1:3) * Sv(1:3,1:3) * B(:, 1:3)';
    
    Ucell = U - gather(neuropil); %what's left over for cell model
    PixL = PixL';
    
%     keyboard;
    
    if k==1
        iclust = initialize_clusters(Ucell, Nk, 'squares', Lx, Ly);       % 'random', Voronoi', 'squares'
    end
    
    % recompute cell activities
    vs = zeros(nSVD, Nk); % initialize cell activities
    for i = 1:Nk
        ix = find(iclust==i);
        if numel(ix)>0
            vs(:,i) = M(ix) * Ucell(:, ix)';
        end
    end
    
    nuVS = sum(vs.^2,1).^.5 + 1e-8;
    vs = bsxfun(@rdivide, vs, nuVS);% normalize activity vectors
    vs = single(vs);
    
    % recompute pixels' assignments
    xs          = vs' * Ucell;
    
    % exclude the pixel's contribution from its own cluster
%             U2 = sum(Ucell.^2,1);
%             xs(iclust + (0:Nk:numel(xs)-1)) = xs(iclust + (0:Nk:numel(xs)-1)) - ...
%                 (M(:) .* U2(:))'./nuVS(iclust);
    %
    [M, iclust] = max(xs,[],1);
    Uneu        = U - bsxfun(@times, M, vs(:,iclust)); %what's left over for neuropil model
    
    %     err(k)      = sum(sum((Uneu-neuropil).^2)).^.5;
    err(k)      = norm(Uneu(:)-neuropil(:));
    
    footPrint = get_footprint(xs, Ly, Lx, ops);
        
    if 1
        %---------------------------------------------%
        % remove clusters
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
        lam = M;
        for i = 1:Nk
            ix = find(iclust==i);
            nT0 = numel(ix);
            if nT0>0
                vM = lam(ix);
                                vM = vM/sum(vM.^2)^.5;
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

keyboard;

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

res.S          = Sm;
res.Sraw       = S;

res.iclust  = iclust;
res.M       = M;
res.lambda  = lam;

%
res.Ly  = Ly;
res.Lx  = Lx;
stat    = get_stat(res);

if isfield(ops, 'ResultsSavePath')
    iplane = ops.iplane;
    if ~exist(ops.ResultsSavePath, 'dir')
        mkdir(ops.ResultsSavePath)
    end

    save(sprintf('%s/F_%s_%s_plane%d_Nk%d.mat', ops.ResultsSavePath, ...
        ops.mouse_name, ops.date, iplane, Nk),  'ops', 'res', 'stat');
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
