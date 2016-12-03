function iclust = initialize_clusters(Ucell, Nk, type, Lx, Ly)

switch type
    case 'random'
        vs = randn(size(Ucell,1), Nk);
        vs = bsxfun(@rdivide, vs, sum(vs.^2,1).^.5 + 1e-8);% normalize activity vectors
        vs = single(vs);
        xs          = vs' * Ucell;
        [~, iclust] = max(xs,[],1);
    case 'Voronoi'
        
        xs = repmat(1:Lx, Ly, 1);
        ys = repmat((1:Ly)', 1, Lx);
        
        randx = rand(1, Nk) * Lx;
        randy = rand(1, Nk) * Ly;
        
        dx = repmat(xs(:), 1, Nk) - repmat(randx, numel(xs(:)), 1);
        dy = repmat(ys(:), 1, Nk) - repmat(randy, numel(ys(:)), 1);
        
        dxy = dx.^2 + dy.^2;
        [~, iclust] = min(dxy, [], 2);
    case 'squares'
        nsqrt = round(sqrt(Nk));
        
        xs = repmat(round(linspace(1, nsqrt, Lx)), Ly, 1);
        ys = repmat(round(linspace(1, nsqrt, Ly))', 1, Lx);
        iclust = xs + (ys-1) * nsqrt;
        
end


iclust = iclust(:)';
