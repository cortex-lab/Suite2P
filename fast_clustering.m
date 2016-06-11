function [ops, stat, res] = fast_clustering(ops, U, Sv)

U =  reshape(U, [], size(U,ndims(U)));
iplane = ops.iplane;

for i = 1:size(U,2)
  U(:,i) = U(:,i)  * Sv(i).^.5;
end
U = U';
[nSVD, Npix] = size(U);
%% clustering options
Ly = numel(ops.yrange);
Lx = numel(ops.xrange);

ops.Nk0 = ceil(sqrt(ops.Nk0)).^2;

Nk      = ops.Nk0;
niter   = ops.niterclustering;



xs = repmat(1:Lx, Ly, 1);
ys = repmat((1:Ly)', 1, Lx);



randx = rand(1, Nk) * Lx;
randy = rand(1, Nk) * Ly;

dx = repmat(xs(:), 1, Nk) - repmat(randx, numel(xs(:)), 1);
dy = repmat(ys(:), 1, Nk) - repmat(randy, numel(ys(:)), 1);

dxy = dx.^2 + dy.^2;

%%
[~, iclust] = min(dxy, [], 2);

clear dx dy

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

M = ones(Npix,1, 'single');
%%
tic
ison = true(Nk,1);

for k = 1:niter  
%     [Nkiter(k) numel(unique(iclust)) sum(ison)]
    
    vs = zeros(nSVD, Nk, 'single');
    for i = 1:Nk
        if ison(i)
            ix = find(iclust==i);
            
            nT0 = numel(ix);
            if nT0>0
                vM = M(ix);
                vs(:,i) = U(:, ix) * vM;
            end
        end
    end
  
    vs = normc(vs);
    xs = vs' * U;
    %
    
    %
    [M, iclust] = max(abs(xs),[],1);    
    xs(iclust + (0:Nk:numel(xs)-1)) = 0;
    [M2, iclust2] = max(abs(xs),[],1);
    
    dM = M - M2;
    indrem = Nkiter(k) - Nkiter(k+1);
    dMk = Inf*ones(Nk,1);
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
    
    
    M       = M';
    iclust  = iclust';
    
    err(k) = sum(M(:));
    if rem(k,10)==1 && ops.ShowCellMap
        %%
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
        V = max(0, min(.5 * reshape(lam, Ly, Lx)/mean(lam(:)), 1));
        H = reshape(r(iclust), Ly, Lx);
        rgb_image = hsv2rgb(cat(3, H, Sat, V));
        imagesc(rgb_image)
        axis off
        drawnow
        fprintf('explained variance is %2.6f time %2.2f \n', err(k), toc)
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
%%
clear res

res.iclust  = iclust;
res.M       = M;
res.lambda  = lam;

%%
res.Ly  = Ly;
res.Lx  = Lx;
stat    = get_stat(res);
%%
if ~exist(ops.ResultsSavePath, 'dir')
    mkdir(ops.ResultsSavePath)
end
save(sprintf('%s/F_%s_%s_plane%d_Nk%d.mat', ops.ResultsSavePath, ...
    ops.mouse_name, ops.date, iplane, Nk),  'ops', 'res', 'stat')

%%
% sk = skewness(F,[],2);
% [~, isk] = sort(sk, 'descend');
% clf
% for i = 1:20
%    plot(5*i + zscore(F(isk(i), :)))
%    hold all
% end
% axis tight
