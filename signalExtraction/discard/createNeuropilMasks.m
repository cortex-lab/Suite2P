% creates surround neuropil masks which exclude the cells
function [ops, neuropMasks]=createNeuropilMasks(ops, stat, cellPix)
%
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


ops.innerNeuropil  = getOr(ops, 'innerNeuropil', 1); % padding around cell to exclude from neuropil
ops.outerNeuropil  = getOr(ops,'outerNeuropil', Inf); % radius of neuropil surround
%

outRadius  = ops.outerNeuropil;
inRadius   = ops.innerNeuropil;

if isinf(outRadius)
    ops.minNeuropilPixels = getOr(ops, 'minNeuropilPixels', 350); 
    ops.ratioNeuropil     = getOr(ops, 'ratioNeuropil', 5);
end


[Ny, Nx]=size(cellPix);

[xx_np, yy_np] = ndgrid((1:Ny),(1:Nx));

%compute the inner border of the ring-like shaped surrounding neuropils
se = strel('disk',inRadius);
expandedGeneralMask = sign(imdilate(cellPix,se));
% expandedGeneralMask = sign(cellPix);

neuropMasks         = zeros(length(stat), Ny, Nx,'single');

%%
%compute the mask of each surrounding neuropils
for k = 1:length(stat)
    centerCell   = [stat(k).med(1) stat(k).med(2)];

    if stat(k).radius > 0
        
        if isinf(ops.outerNeuropil)
            radiusCell         = stat(k).radius;
            outRadius          = ops.ratioNeuropil * radiusCell;
            outRadius          = max(outRadius, sqrt(ops.minNeuropilPixels/pi));
            totPix             = 0;
            its                = 0;
            while totPix < ops.minNeuropilPixels
                neuropRegion       = sqrt((xx_np-centerCell(1)).^2 + ...
                    (yy_np-centerCell(2)).^2)<=outRadius;
                neuropNoCells      = neuropRegion - expandedGeneralMask > 0;
                totPix             = sum(neuropNoCells(:));
                outRadius          = outRadius * 1.25;
            end
            
            % embed neuropil masks in full FOV (LyU x LxU)
            neuropMasks(k,:,:) = single(neuropNoCells)/totPix;
        else
            neuropRegion       = sqrt((xx_np-centerCell(1)).^2 + ...
                (yy_np-centerCell(2)).^2)<=outRadius;
            neuropNoCells      = neuropRegion - expandedGeneralMask > 0;
            totPix             = sum(neuropNoCells(:));
            neuropMasks(k,:,:) = single(neuropNoCells)/totPix;
        end
        
    end
end

