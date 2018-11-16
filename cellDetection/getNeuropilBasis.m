% compute neuropil basis functions for cell detection (either cosyne or fourier)
function [S,ops] = getNeuropilBasis(ops, Ly, Lx, type)

ops.ratioNeuropil = getOr(ops, 'ratioNeuropil', 3);


nTilesY = ceil(Ly / (ops.ratioNeuropil * ops.diameter/2)); % neuropil is modelled as nTiles by nTiles
nTilesX = ceil(Lx / (ops.ratioNeuropil * ops.diameter/2)); % neuropil is modelled as nTiles by nTiles
% for ops.diameter = 8 and 512x512 FOV, default is 22 x 22 tiles
S = zeros(Ly, Lx, nTilesY, nTilesX, 'single');

xs = 1:Lx;
ys = 1:Ly;
        
switch type
    case 'raisedcosyne'
        xc = linspace(1, Lx, nTilesX);
        yc = linspace(1, Ly, nTilesY);
        yc = yc';
        
        sigx = 4*(Lx - 1)/nTilesX;
        sigy = 4*(Ly - 1)/nTilesY;
        
        
        for kx = 1:nTilesX
            for ky = 1:nTilesY
                cosx = 1+cos(2*pi*(xs - xc(kx))/sigx);
                cosy = 1+cos(2*pi*(ys - yc(ky))/sigy);
                cosx(abs(xs-xc(kx))>sigx/2) = 0;
                cosy(abs(ys-yc(ky))>sigy/2) = 0;                                
                S(:, :,ky, kx) = cosy' * cosx;
            end
        end
    case 'dct'
        Ay = dctbasis(Ly,1);
        Ax = dctbasis(Lx,1);
        
        for j = 1:nTilesX
            for i = 1:nTilesY
                S(:,:,i,j) = Ay(:,j) * Ax(:,i)';
            end
        end
    case 'Fourier'
        nTilesY = 1 + 2*ceil(nTilesY/2);
        nTilesX = 1 + 2*ceil(nTilesX/2);
        S = zeros(Ly, Lx, nTilesY, nTilesX, 'single');
        
        Ay = ones(Ly, nTilesY);
        Ax = ones(Lx, nTilesX);
        for k = 1:(nTilesY-1)/2
            Ay(:,2*k)   = sin(2*pi*ys*k/Ly);
            Ay(:,2*k+1) = cos(2*pi*ys*k/Ly);
        end
        for k = 1:(nTilesX-1)/2
            Ax(:,2*k)   = sin(2*pi*xs*k/Lx);
            Ax(:,2*k+1) = cos(2*pi*xs*k/Lx);
        end
        
        for j = 1:nTilesX
            for i = 1:nTilesY
                S(:,:,i,j) = Ay(:,j) * Ax(:,i)';
            end
        end
end

S = reshape(S, [], nTilesX*nTilesY);
S = normc(S);