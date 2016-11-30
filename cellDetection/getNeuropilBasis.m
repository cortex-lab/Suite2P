function S = getNeuropilBasis(ops, type)

Ly = ops.Ly;
Lx = ops.Lx;

switch type
    case 'raisedcosyne'
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
end