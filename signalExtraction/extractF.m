
% computes cell and neuropil fluorescence for surround model of neuropil
function [ops, stat, Fcell, FcellNeu] = extractF(ops, stat)

Nk       = numel(stat); % all ROIs

[Ly, Lx] = size(ops.mimg);
% stat     = getNonOverlapROIs(stat, Ly, Lx);

% create cell masks and cell exclusion areas
[stat, cellPix, cellMasks] = createCellMasks(stat, Ly, Lx);

% create surround neuropil masks
[ops, neuropMasks] = neuropilMasks2(ops, stat, cellPix);

% add surround neuropil masks to stat
for k = 1:Nk
    stat(k).ipix_neuropil = find(squeeze(neuropMasks(k,:,:))>0);
end

% convert masks to sparse matrices for fast multiplication
% neuropMasks = sparse(double(neuropMasks(:,:)));
% cellMasks   = sparse(double(cellMasks(:,:)));
neuropMasks = single(neuropMasks(:,:));
cellMasks   = single(cellMasks(:,:));
% get fluorescence and surround neuropil
nimgbatch = 1000;
ix = 0;
fclose all;
fid = fopen(ops.RegFile, 'r');

tic
F = zeros(Nk, sum(ops.Nframes), 'single');
Fneu = zeros(Nk, sum(ops.Nframes), 'single');

while 1
    data = fread(fid,  Ly*Lx*nimgbatch, '*int16');
    if isempty(data)
        break;
    end
    
    data = reshape(data, Ly, Lx, []);    
    NT   = size(data,3);
    data = reshape(data, [], NT);
    data = single(data);
    
    % process the data
    %data = my_conv2(data, ops.sig, [1 2]);
    
    % compute cell fluorescence
    % each mask is weighted by lam (SUM TO 1)
    F(:,ix + (1:NT)) = cellMasks * data;
    
    % compute neuropil fluorescence
    Fneu(:,ix + (1:NT)) = neuropMasks * data;
    
    ix = ix + NT;
    if rem(ix, 3*NT)==0
        fprintf('Frame %d done in time %2.2f \n', ix, toc)
    end
end
fclose(fid);

% keyboard;

% get activity stats
[stat, F, Fneu] = getActivityStats(ops, stat, F, Fneu);

%
csumNframes = [0 cumsum(ops.Nframes)];
Fcell       = cell(1, length(ops.Nframes));
FcellNeu    = cell(1, length(ops.Nframes));
for i = 1:length(ops.Nframes)
    Fcell{i}     = F(:, csumNframes(i) + (1:ops.Nframes(i)));
    FcellNeu{i}  = Fneu(:, csumNframes(i) + (1:ops.Nframes(i)));
end

end

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
        if sum(cell_pix(ypix1 + (xpix1-1)*Ly)<.5) - nring > ops.minNeuropilPixels
            break;
        end
    end    
    ix = cell_pix(ypix1 + (xpix1-1)*Ly)<.5;
    ypix1 = ypix1(ix);
    xpix1 = xpix1(ix);
    neuropMasks(k,ypix1+(xpix1-1)*Ly) = 1;
    neuropMasks(k,ypix+(xpix-1)*Ly)   = 0;
end
neuropMasks = neuropMasks ./ sum(sum(neuropMasks,2),3);
end

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
        stat(k).radius = max(3, stat(k).radius);
        
        % put cell mask into cellMasks for computing fluorescence
        cellMasks(k, ipix) = lam / sum(lam);
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

cellPix = min(1, cellPix);
end

function [ypix, xpix] = extendROI(ypix, xpix, Ly, Lx, niter)
for k = 1:niter
    yx = [ypix xpix; ypix xpix+1; ypix xpix-1; ypix-1 xpix; ypix+1 xpix];
    yu = unique(yx, 'rows');
    ix = yu(:,1)>0 & yu(:,1)<=Ly &  yu(:,2)>0 & yu(:,2)<=Lx;
    ypix = yu(ix,1);
    xpix = yu(ix,2);
end
end
