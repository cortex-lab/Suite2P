% creates cell masks for fluorescence computation
% also creates exclusion mask for neuropil (cellPix)
function [stat, cellPix, cellMasks] = createCellMasks(stat, Ny, Nx, allow_overlap)

if nargin > 3
    aov = allow_overlap;
else
    aov = 0;
end

Nk = length(stat);
cellPix = zeros(Ny, Nx);
cellMasks = zeros(Nk, Ny, Nx, 'single');
for k = 1:Nk
    % only use non-overlapping pixels of cell
<<<<<<< HEAD
    ipix = stat(k).ipix(stat(k).isoverlap==0);
    ypix = stat(k).ypix(stat(k).isoverlap==0);
    xpix = stat(k).xpix(stat(k).isoverlap==0);
    lam  = stat(k).lam(stat(k).isoverlap==0);
    
    stat(k).radius = NaN;
    if ~isempty(ipix)
        % fit MV gaussian to cell mask
        % define cell as all pixels within 2 std's of lambda
        params      = FitMVGaus(ypix, xpix, lam, 2);
        
        % cell radius
        stat(k).radius = sqrt(mean(params.eval));
        
        % put cell mask into cellMasks for computing fluorescence
        cellMasks(k, ipix) = lam / sum(lam);
    end
    % use all pixels for neuropil mask computation
    temp        = zeros(Ny, Nx);
    ipix = stat(k).ipix;
    lam  = stat(k).lam;
    temp(ipix)  = lam;
    
    % input thresholded pixels to cellPix to exclude in neuropil computation
    cellPix = cellPix + (temp > 0);
end

=======
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
    
    temp(ipix)  = lam;
    
    % fit MV gaussian to cell mask
    % define cell as all pixels within 2 std's of lambda
    if ~isempty(ypix)
        params      = FitMVGaus(ypix, xpix, lam, 2);
        % cell radius
        stat(k).radius = sqrt(mean(params.eval));
        xy          = max(1, min(repmat([Ny Nx],size(params.xy,1),1), round(params.xy)));
        ithres      = sub2ind([Ny Nx], xy(:,1), xy(:,2));
        lamthres    = mean(temp(ithres));
        temp        = zeros(Ny, Nx);
        temp(ipix(lam > lamthres))  = lam(lam > lamthres);
    
        % input thresholded pixels to cellPix to exclude in neuropil computation
        cellPix = cellPix + (temp > 0);
    
        % put cell mask into cellMasks for computing fluorescence
        cellMasks(k, ipix) = lam / sum(lam);
    else
        stat(k).radius = 0;
    end
end

cellPix = min(1, cellPix);
>>>>>>> 3d5bdb6d52d4602bcf83c766950efeace8bcb498
