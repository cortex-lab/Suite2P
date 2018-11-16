% creates surround neuropil masks which exclude the cells
function [ops, neuropMasks] = neuropilMasks2(ops, stat, cell_pix)

% THE DESCRIPTION BELOW IS NO LONGER QUITE RIGHT
% ops.innerNeuropil: distance between cell mask and inner boundary of neuropil
% mask (in pixels), default is 3 pixels
% IF NEUROPIL FUNCTION OF CELL SIZE
% ops.minNeuropilPixels: minimum number of pixels in neuropil surround
% ops.ratioNeuropil: radius of surround / radius of cell
% IF NEUROPIL FIXED SIZE
% ops.outerNeuropil: outer radius of Neuropil mask -- if set to Inf then
% neuropil is a function of the cell size
%
% OUTPUT:
% neuropMasks: [Nk x Ny x Nx] matrix identify different neuropil masks with an
% integer index
%

[Ly, Lx] = size(cell_pix);
ops.innerNeuropil  = getOr(ops, 'innerNeuropil', 2); % padding around cell to exclude from neuropil

neuropMasks         = zeros(length(stat), Ly, Lx,'single');

%compute the mask of each surrounding neuropils
for k = 1:length(stat)
    ypix = stat(k).ypix;
    xpix = stat(k).xpix;
    % first extend to get ring of dis-allowed pixels
    [ypix, xpix] = extendROI(ypix, xpix, Ly, Lx, ops.innerNeuropil);
    nring = sum(cell_pix(ypix + (xpix-1)*Ly)<.5);
    ypix1 = ypix;
    xpix1 = xpix;
    for j = 1:100
        [ypix1, xpix1] = extendROI(ypix1, xpix1, Ly, Lx, 5);
        if sum(cell_pix(ypix + (xpix-1)*Ly)<.5) - nring > ops.minNeuropilPixels
            break;
        end
    end    
    ix = cell_pix(ypix + (xpix-1)*Ly)<.5;
    ypix1 = ypix1(ix);
    xpix1 = xpix1(ix);
    neuropMasks(k,ypix1+(xpix1-1)*Ly) = 1;
    neuropMasks(k,ypix+(xpix-1)*Ly)   = 0;
end

neuropMasks = neuropMasks ./ sum(sum(neuropMasks,2),3);



end