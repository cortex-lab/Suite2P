% run cell detection on spatial masks U
function [ops, stat] = sourcery_inc(ops)

tic
[ops, U0, U2]    = get_svdForROI(ops);
        

ops.fig         = getOr(ops, 'fig', 1);
ops.ThScaling   = getOr(ops, 'ThScaling', 1);

% reshape U to be (nMaps x Y x X)
U0 = permute(U0, [3 1 2]);
[nSVD, Ly, Lx] = size(U0);
% [nSVD, Npix] = size(U0);
% U0 = reshape(U0, nSVD, Ly, Lx);

% compute neuropil basis functions for cell detection
[S, ops] = getNeuropilBasis(ops, Ly, Lx, 'Fourier'); % 'raisedcosyne', 'Fourier'
S = normc(S)';
nBasis = size(S,1);

StU = U0(:,:) * S'; % covariance of neuropil with spatial masks
StS = S*S'; % covariance of neuropil basis functions

% make cell mask with ops.diameter
d0   = ops.diameter; % expected cell diameter in pixels
sig = ceil(ops.diameter/4); 

iter = 0;
icell = 0;
r = rand(1e4,1);
%%
L   = zeros(0, Ly*Lx, 'single');
LtU = zeros(nSVD, 0, 'single');
LtS = zeros(nBasis, 0, 'single');

% regress maps onto basis functions and subtract neuropil contribution
% U = Ucell + neu'*S'
% neu = inv(S'*S) * (S'*U')
neu     = StS\StU';
Ucell   =  U0 - reshape(neu' * S, size(U0));

refine = -1;

%
tic
while 1
    iter = iter + 1;   
    if refine<0
        % residual is smoothed at every iteration
        us = my_conv2_circ(Ucell, sig, [2 3]);
        % compute log variance at each location
        V = sq(mean(us.^2,1));
        V = double(V);
        um = sq(mean(Ucell.^2,1));
        um = my_conv2_circ(um, sig, [1 2]);
        V = V./um ;
        %     V = log(V./um);
        V = double(V);
        % do the morphological opening trick
        % take the running max of the running min
        % this normalizes the brightness of the image
        if iter==1
            lbound = -my_min2(-my_min2(V, d0), d0);
        end
        V = V - lbound;
        
        if iter==1
            % find indices of all maxima  in plus minus 1 range
            % use the median of these peaks to decide stopping criterion
            maxV    = -my_min(-V, 1, [1 2]);
            ix      = (V > maxV-1e-10);
            % threshold is the mean peak, times a potential scaling factor
            pks = V(ix);
            Th  = ops.ThScaling * median(pks(pks>1e-4));
            ops.Vcorr = V;
        end
        
        % just in case this goes above original value
        V = min(V, ops.Vcorr);
        
        ncells = icell;
        
        
        while icell<ncells+200
            [~, ind] = max(V(:));
            [i, j] = ind2sub(size(V), ind);
            if V(i,j) < Th
                break;
            end
            
            icell = icell + 1;
            [yp, xp, la, ix, code] = iter_extend(i,j, Ucell, us(:, i,j), -1);
            ypix{icell} = yp;
            xpix{icell} = xp;
            lam{icell} = la;
            
            Ucell(:, yp+(xp-1)*Ly) = Ucell(:, yp+(xp-1)*Ly) - code * la';
            
            [yp, xp] = extendROI(yp, xp, Ly, Lx, round(mean(d0)));
            V(yp + (xp-1)*Ly) = 0;
        end
        
        newcells = icell-ncells;
        if iter==1
            Nfirst = newcells;
        end        
        
        LtU(nSVD, icell) = 0;
        LtS(nBasis, icell) = 0;
        L(icell, 1) = 0;
        for n = ncells+1:icell
            L(n, ypix{n} + (xpix{n}-1)*Ly) = lam{n};
            LtU(:, n) = U0(:, ypix{n}+(xpix{n}-1)*Ly) * lam{n};
            LtS(:, n) = S(:, ypix{n}+(xpix{n}-1)*Ly) * lam{n};
        end
        
        % ADD NEUROPIL INTO REGRESSION HERE
        LtL     = full(L*L');
        codes   = ([LtL LtS'; LtS StS]+ 1e-3 * eye(icell+nBasis))\[LtU StU]';
        neu     = codes(icell+1:end,:);
        codes   = codes(1:icell,:);
        %     codes = (LtL+ 1e-3 * eye(icell))\LtU;
        
        ncells = icell;
    end
    
    % subtract off everything
    Ucell = U0 - reshape(neu' * S, size(U0)) - reshape(double(codes') * L, size(U0));    
    
    % re-estimate masks
    n = 1;    
    for k = 1:ncells        
        Ucell(:, ypix{n}+(xpix{n}-1)*Ly) = Ucell(:, ypix{n}+(xpix{n}-1)*Ly) + codes(k,:)' * lam{n}';
        [ypix{n}, xpix{n}, lam{n}, ix, codes(n,:)] = iter_extend(ypix{n}, xpix{n}, Ucell, codes(k,:)', refine);
        if sum(ix)==0
           disp('dropped ROI with no pixels') 
           ypix(n)= [];
           xpix(n)= [];
           lam(n)= [];
           continue;
        end        
        Ucell(:, ypix{n}+(xpix{n}-1)*Ly) = Ucell(:, ypix{n}+(xpix{n}-1)*Ly) - codes(n,:)' * lam{n}';        
        n = n+1;
    end
    n = n-1;
    codes = codes(1:n, :);
    ncells = n;
    
    L = zeros(ncells, Ly*Lx, 'single');
    LtU = zeros(nSVD, ncells, 'single');
    LtS = zeros(nBasis, ncells, 'single');    
    for n = 1:ncells
        L(n, ypix{n} + (xpix{n}-1)*Ly) = lam{n};
        if refine<0
            LtU(:, n) = U0(:, ypix{n}+(xpix{n}-1)*Ly) * lam{n};
            LtS(:, n) = S(:, ypix{n}+(xpix{n}-1)*Ly) * lam{n};
        end
    end
    
    err(iter) = mean(Ucell(:).^2);
    fprintf('time %2.2f, %d total ROIs, err %4.4f, thresh %4.4f \n', toc, icell, err(iter), Th)
    if ops.fig           
        figure(13132)
        subplot(1,2, 1);
        imagesc(ops.Vcorr, [0 2*Th])
        axis off        
        subplot(1,2, 2);
        imagesc(V, [0 2*Th])
        axis off        
        figure(13133)
        stat = struct('ypix', ypix, 'xpix', xpix, 'lam', lam);
        drawClusters(ops, stat, r)        ;
        drawnow
    end
    
    if refine==0
        break;
    end
    if refine==2
        stat = struct('ypix', ypix, 'xpix', xpix, 'lam', lam);
        stat = connectedRegion2(stat, ops);
        [stat, ix] = removeOverlaps(stat, ops, Ly, Lx);
        fprintf('removed %d overlapping ROIs\n', numel(ypix)-numel(ix));        
        ypix = {stat.ypix};
        xpix = {stat.xpix};
        lam = {stat.lam};
        L = L(ix, :, :);
        codes = codes(ix, :);
        ncells = numel(ypix);    
    end
    if refine>0
        Ucell = Ucell + reshape(neu' * S, size(U0));    
    end
    if refine<0 && (newcells<Nfirst * getOr(ops, 'stopSourcery', 1/10)) || (iter>= getOr(ops, 'maxIterRoiDetection', 20))
        refine = 3;
        U0 = permute(U2, [3 1 2]);
        Ucell = U0;
    end
    if refine>=0
        StU = Ucell(:,:) * S'; % covariance of neuropil with spatial masks
        neu     = StS\StU';
    end
    refine = refine - 1;    
end

% subtract off neuropil only
Ucell = U0 - reshape(neu' * S, size(U0));

stat = struct('ypix', ypix, 'xpix', xpix, 'lam', lam);
for j = 1:length(stat)
   stat(j).lam = stat(j).lam .* ops.sdmov(stat(j).ypix + (stat(j).xpix-1)*Ly);
end

stat = postprocess(ops, stat, Ucell, codes);

% stat = getFootprint(ops, codes, Ucell, mPix, mLam);
% stat = anatomize(ops, mPix, mLam, stat);
% stat = weightsMeanImage(ops, stat, model);

end

function [stat, ix] = removeOverlaps(stat, ops, Ly, Lx)
ops.max_overlap = getOr(ops, 'max_overlap', 0.75);
ncells = numel(stat);
mask = zeros(Ly,Lx);
ix = 1:ncells;
for n = 1:ncells
   ypix = stat(n).ypix; 
   xpix = stat(n).xpix;
   mask(ypix + (xpix-1)*Ly) = mask(ypix + (xpix-1)*Ly) + 1;
end
while 1
    O = zeros(numel(stat), 1);
    for n = 1:numel(stat)
        ypix = stat(n).ypix;
        xpix = stat(n).xpix;
        O(n) = mean(mask(ypix + (xpix-1)*Ly) > 1.5);
    end
    [~, i] = max(O);
    if O(i) > ops.max_overlap
        ypix = stat(i).ypix;
        xpix = stat(i).xpix;
        mask(ypix + (xpix-1)*Ly)  = mask(ypix + (xpix-1)*Ly) - 1;
        stat(i) = [];
        ix(i) = [];
    else
        break;
    end
end
end

function stat = getConnected(Ly, Lx, stat)
ypix = stat.ypix;
xpix = stat.xpix;
lam = stat.lam;
[~, i0] = max(lam);
mask = zeros(Ly,Lx);
mask(ypix + (xpix-1)*Ly) = lam;
ypix = ypix(i0);
xpix = xpix(i0);
nsel = 1;
while 1
   [ypix, xpix] = extendROI(ypix, xpix, Ly, Lx, 1) ;
   ix = mask(ypix + (xpix-1)*Ly) > 1e-10;
   ypix = ypix(ix);
   xpix = xpix(ix);
   if numel(ypix) <= nsel
       break;
   end
   nsel = numel(ypix);
end
lam = mask(ypix + (xpix-1)*Ly);
stat.ypix = ypix;
stat.xpix = xpix;
stat.lam = lam;
end

function stat = getOverlaps(stat, ops)
Ly = ops.Ly; 
Lx = ops.Lx; 
ncells = numel(stat);
mask = zeros(Ly, Lx);
for n = 1:ncells
    ypix = stat(n).ypix;
    xpix = stat(n).xpix;
    mask(ypix + (xpix-1)*Ly) = mask(ypix + (xpix-1)*Ly) + 1;
end
for n = 1:ncells
    ypix = stat(n).ypix;
    xpix = stat(n).xpix;
    stat(n).isoverlap = mask(ypix + (xpix-1)*Ly) > 1.5;    
end
end

function stat = connectedRegion2(stat, ops)
connected = getOr(ops, 'connected', 1);
Ly = numel(ops.yrange);
Lx = numel(ops.xrange);
if connected
   for j = 1:numel(stat)
      stat(j) = getConnected(Ly, Lx, stat(j)); 
   end
end
end

function stat = postprocess(ops, stat, Ucell, codes)
stat = connectedRegion2(stat, ops); % connected region
stat = getStat(ops, stat, Ucell, codes);
stat = getOverlaps(stat,ops); % remove overlaps
end


function stat = getStat(ops, stat, Ucell, codes)
d0   = ceil(ops.diameter); % expected cell diameter
dx = repmat([-2*d0:2*d0], 4*d0+1, 1);
dy = dx';

rs = (dx.^2 + dy.^2).^.5;
rs = rs(rs<=2*d0);
r2sort = sort(rs);

frac = [.5]; %[0.15 0.25 0.33 0.5 .75];
Ny = length(ops.yrange);
Nx = length(ops.xrange);
Ly = ops.Ly;

% find maximum contamination distance for each ROI
for j = 1:numel(stat)    
    stat(j).lambda  = stat(j).lam;    
    stat(j).med     = [median(stat(j).ypix) median(stat(j).xpix)];
    stat(j).neuropilCoefficient = 0.7;
    stat(j).baseline            = 0;        
    
    % compute "footprint" length of ROI
    y0 = ceil(stat(j).med(1));
    x0 = ceil(stat(j).med(2));    
    
    % pixels within FOV and within radius d0 of center
    [yp, xp] = extendROI(stat(j).ypix, stat(j).xpix, Ny, Nx, d0);    
    rs0 = ((yp-y0).^2 + (xp-x0).^2).^.5;
    ipix = yp + Ny * (xp-1);    
    % weight of cell onto each pixel across components
    proj = codes(j,:) * Ucell(:, ipix); 
    % mean of radius of pixels with weights > frac
    stat(j).footprint = mean(rs0(proj>max(proj) * frac)) /ops.diameter;    
    
    params = FitMVGaus(stat(j).ypix,stat(j).xpix, stat(j).lam);
    stat(j).aspect_ratio = sqrt(max(params.eval) / min(params.eval));
    stat(j).ellipse      = params.xy;    
    
    r2 = (stat(j).ypix-y0).^2 + (stat(j).xpix-x0).^2;
    r2 = sqrt(r2);
    stat(j).mrs  = mean(r2);
    stat(j).mrs0  = mean(r2sort(1:min(numel(r2sort), numel(r2))));    
    stat(j).compact = stat(j).mrs/stat(j).mrs0;
    
    stat(j).ypix = stat(j).ypix + ops.yrange(1)-1;
    stat(j).xpix = stat(j).xpix + ops.xrange(1)-1;
    stat(j).ipix    = stat(j).ypix + Ly * (stat(j).xpix-1);
    stat(j).npix    = numel(stat(j).ipix);
    stat(j).med     = [median(stat(j).ypix) median(stat(j).xpix)];
end

mpix = mean([stat.npix]);
mfoot = median([stat.footprint]);
for j = 1:numel(stat)
    stat(j).footprint = stat(j).footprint/mfoot;
    stat(j).npix_norm = stat(j).    npix/mpix;
end
end

function [ypix, xpix, lam, ix, code] = iter_extend(ypix, xpix, Ucell, code, refine)
[nsvd, Lyc, Lxc] = size(Ucell);
npix = 0;

while npix<10000    
    npix = numel(ypix);
    [ypix, xpix] = extendROI(ypix,xpix,Lyc,Lxc, 1);
    usub = Ucell(:, ypix+(xpix-1)*Lyc);
    lam = code' * usub;
    lam = lam(:); % squeeze only index 2
    ix = lam > max(0, max(lam)/5);
    if sum(ix)==0
        break;
    end
    ypix = ypix(ix);
    xpix = xpix(ix);
    lam = lam(ix);
    lam = lam/sum(lam.^2+1e-10)^.5;
    if refine<0
        code = usub(:, ix) * lam;
    end
    if npix>=sum(ix)
        break;
    else
        npix = numel(ypix);
    end    
end
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

function drawClusters(ops, stat, r)
Ly = numel(ops.yrange);
Lx = numel(ops.xrange);
ncells = length(stat);
r = r(1:ncells);
Lam = zeros(Ly,Lx);
H = zeros(Ly,Lx);
for n = 1:ncells
    ypix = stat(n).ypix;
    xpix = stat(n).xpix;
    lam  = stat(n).lam;
    isingle = Lam(ypix + (xpix-1)*Ly) + 1e-4 < lam;
    y = ypix(isingle);
    x = xpix(isingle);
    Lam(y + (x-1)*Ly) = lam(isingle);
    H(y + (x-1)*Ly) = r(n) * ones(numel(y), 1);
end

Sat = ones(Ly,Lx);
V = max(0, min(1, .75 * Lam / mean(Lam(Lam>1e-10))));
rgb_image = hsv2rgb(cat(3, H, Sat, V));
imagesc(rgb_image)
axis off
end
