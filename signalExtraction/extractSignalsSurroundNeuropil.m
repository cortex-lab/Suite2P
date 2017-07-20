% computes cell and neuropil fluorescence for surround model of neuropil
function [ops, stat, Fcell, FcellNeu] = extractSignalsSurroundNeuropil(ops, stat)

% IF NEUROPIL FUNCTION OF CELL SIZE
% ops.minNeuropilPixels: minimum number of pixels in neuropil surround
% ops.ratioNeuropil: radius of surround / radius of cell
% IF NEUROPIL FIXED SIZE
% ops.outerNeuropil: outer radius of Neuropil mask -- if set to Inf then
% neuropil is a function of the cell size
ops.innerNeuropil = getOr(ops, 'innerNeuropil', 1);
ops.outerNeuropil = getOr(ops, 'outerNeuropil', Inf);
if isinf(ops.outerNeuropil)
    ops.minNeuropilPixels = getOr(ops, 'minNeuropilPixels', 400);
    ops.ratioNeuropil     = getOr(ops, 'ratioNeuropil', 5);
end

Nk       = numel(stat); % all ROIs

Ny = numel(ops.yrange);
Nx = numel(ops.xrange);

stat = getNonOverlapROIs(stat, Ny, Nx);

[Ly, Lx] = size(ops.mimg);

% create cell masks and cell exclusion areas
[stat, cellPix, cellMasks] = createCellMasks(stat, Ny, Nx);

% create surround neuropil masks
neuropMasks=createNeuropilMasks(ops, stat, cellPix);

% add surround neuropil masks to stat
for k = 1:Nk
    stat(k).ipix_neuropil = find(squeeze(neuropMasks(k,:,:))>0);
end

% convert masks to sparse matrices for fast multiplication
neuropMasks = sparse(double(neuropMasks(:,:)));
cellMasks   = sparse(double(cellMasks(:,:)));

%% get fluorescence and surround neuropil
nimgbatch = 2000;
ix = 0;
fclose all;
fid = fopen(ops.RegFile, 'r');

tic
F = NaN(Nk, sum(ops.Nframes), 'single');
Fneu = NaN(Nk, sum(ops.Nframes), 'single');
while 1
    data = fread(fid,  Ly*Lx*nimgbatch, '*int16');
    if isempty(data)
        break;
    end
    
    %
    data = reshape(data, Ly, Lx, []);
    data = data(ops.yrange, ops.xrange, :);
    data = double(data);
    NT   = size(data,3);
    data = reshape(data, [], NT);
    
    % process the data
    %data = my_conv2(data, ops.sig, [1 2]);
    
    % compute cell fluorescence
    % each mask is weighted by lam (mean 1)
    Ftemp            = cellMasks * data;
    F(:,ix + (1:NT)) = Ftemp;
    
    % compute neuropil fluorescence
    Ftemp               = neuropMasks * data;
    Fneu(:,ix + (1:NT)) = Ftemp;
    
    ix = ix + NT;
    if rem(ix, 3*NT)==0
        fprintf('Frame %d done in time %2.2f \n', ix, toc)
    end
end
fclose(fid);

%%
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
