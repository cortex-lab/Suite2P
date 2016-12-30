function S = getNeuropilBasis(ops, Ly, Lx, type)

TileFactor = getOr(ops, {'TileFactor'}, 1); % this option can be overwritten by the user
nTiles = ceil(TileFactor * (Ly+Lx)/2 / (3 * ops.diameter)); % neuropil is modelled as nTiles by nTiles

S = zeros(Ly, Lx, nTiles, nTiles, 'single');

switch type
    case 'raisedcosyne'
        xc = linspace(1, Lx, nTiles);
        yc = linspace(1, Ly, nTiles);
        yc = yc';
        xs = 1:Lx;
        ys = 1:Ly;
        
        sigx = 4*(Lx - 1)/nTiles;
        sigy = 4*(Ly - 1)/nTiles;
        
        
        for kx = 1:nTiles
            for ky = 1:nTiles
                cosx = 1+cos(2*pi*(xs - xc(kx))/sigx);
                cosy = 1+cos(2*pi*(ys - yc(ky))/sigy);
                cosx(abs(xs-xc(kx))>sigx/2) = 0;
                cosy(abs(ys-yc(ky))>sigy/2) = 0;
                
                S(:, :,ky, kx) = cosy' * cosx;
            end
        end
    case 'Fourier'
        Ay = dctbasis(Ly,1);
        Ax = dctbasis(Lx,1);
        
        for j = 1:nTiles
            for i = 1:nTiles
                S(:,:,i,j) = Ay(:,j) * Ax(:,i)';
            end
        end
end

S = reshape(S, [], nTiles^2);
S = normc(S);