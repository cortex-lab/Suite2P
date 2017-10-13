% creates cell masks for fluorescence computation
% also creates exclusion mask for neuropil (cellPix)
function [stat, centerMasks, surroundMasks, neuropMasks] = ...
    createCenterSurroundMasks2(ops, stat, Ny, Nx, allow_overlap, radius0)

if nargin > 4
    aov = allow_overlap;
else
    aov = 0;
end

if nargin > 4
    radius = radius0;
else
    radius = 3;
end

Nk = length(stat);
cellPix = zeros(Ny, Nx);
centerMasks = zeros(Nk, Ny, Nx, 'single');
surroundMasks = zeros(Nk, Ny, Nx, 'single');

for k = 1:Nk
    % only use non-overlapping pixels of cell
    temp = zeros(Ny, Nx);
    if aov
        ipix = stat(k).ipix(stat(k).isoverlap==0);
        ypix = stat(k).ypix(stat(k).isoverlap==0);
        xpix = stat(k).xpix(stat(k).isoverlap==0);
        lam  = stat(k).lam(stat(k).isoverlap==0);
    else
        ipix = stat(k).ipix;
        ypix = stat(k).ypix;
        xpix = stat(k).xpix;
        lam  = stat(k).lam;
    end
        
    % fit MV gaussian to cell mask
    % define cell as all pixels within 2 std's of lambda
    if ~isempty(ypix)
        params      = FitMVGaus(ypix, xpix, lam, 2);
        % cell radius
        stat(k).radius = sqrt(mean(params.eval));
        
        % add masks for cells with radius greater than
        temp        = zeros(Ny, Nx);
        tempS       = zeros(Ny, Nx);
        if stat(k).radius >= radius
            pdist = ((stat(k).ypix - stat(k).med(1)).^2 + ...
                (stat(k).xpix - stat(k).med(2)).^2);
            temp(stat(k).ipix( pdist   < (stat(k).radius*.7)^2)) = 1;
            tempS(stat(k).ipix(pdist >= (stat(k).radius*.7)^2)) = 1;
        end
        
        centerMasks(k, :, :)      = temp / sum(temp(:));
        surroundMasks(k, :, :)    = tempS / sum(tempS(:));
    else
        stat(k).radius = 0;
    end
    
    % use all pixels for neuropil mask computation
    temp        = zeros(Ny, Nx);
    ipix = stat(k).ipix;
    lam  = stat(k).lam;
    temp(ipix)  = lam;
    
    % input thresholded pixels to cellPix to exclude in neuropil computation
    cellPix = cellPix + (temp > 0);
end

% create surround neuropil masks
[~, neuropMasks] = createNeuropilMasks(ops, stat, cellPix);
