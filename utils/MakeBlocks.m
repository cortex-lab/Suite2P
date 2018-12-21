% splits FOV in X and Y as specified by numBlocks and blockFrac
% computes xyMask with smoothing constants specified to smooth registration
% offsets across blocks
function ops = MakeBlocks(ops)

Lx = ops.Lx;
Ly = ops.Ly;
numBlocks = ops.numBlocks;

% fraction of block of total image
bfrac     = 1./max(2,(numBlocks-3));
bfrac(numBlocks==1) = 1;
ops.blockFrac = getOr(ops, {'blockFrac'}, bfrac);
bpix      = round(ops.blockFrac .* [Ly Lx]);

% pixel overlap of blocks
ops.pixoverlap    = [];
ops.pixoverlap = round((bpix.*numBlocks-[Ly Lx])./(numBlocks-1.9));

% edges of blocks
yB        = linspace(0, Ly, numBlocks(1)+1);
yB        = round((yB(1:end-1) + yB(2:end)) / 2);

xB        = linspace(0, Lx, numBlocks(2)+1);
xB        = round((xB(1:end-1) + xB(2:end)) / 2);

% compute indices of blocks
ib        = 0;
for iy = 1:numBlocks(1)
    if iy == numBlocks(1)
        yB(iy)  = Ly - floor(bpix(1)/2);
    elseif iy == 1
        yB(iy)  = floor(bpix(1)/2) + 1;
    end
    if numBlocks(2) > 1
        for ix = 1:numBlocks(2)
            ib = ib+1;
            if ix == numBlocks(2)
                xB(ix)  = Lx - floor(bpix(2)/2);
            elseif ix == 1
                xB(ix)  = floor(bpix(2)/2) + 1;
            end
        
            ops.yBL{ib} = [max(1,yB(iy)-floor(bpix(1)/2)) : ...
                min(Ly,yB(iy)+floor(bpix(1)/2))];
            ops.xBL{ib} = [max(1,xB(ix)-floor(bpix(2)/2)) : ...
                min(Lx,xB(ix)+floor(bpix(2)/2))];
        end
    else
        ib = ib+1;
        ops.yBL{ib} = [max(1,yB(iy)-floor(bpix(1)/2)) : ...
            min(Ly,yB(iy)+floor(bpix(1)/2))];
        ops.xBL{ib} = [1:Lx];
    end
end
nblocks = ib;

% compute smoothing mask xyMask
% gaussian centered on block with width 2/3 the size of block
% (or specified by user as ops.smoothBlocks)
sT(1)        = mean(diff(yB)) * 2/3;
sT(2)        = mean(diff(xB)) * 2/3;
ops.smoothBlocks = getOr(ops, {'smoothBlocks'}, sT);
ops.smoothBlocks = max(10, ops.smoothBlocks);
sT           = ops.smoothBlocks;

xyMask = zeros(Ly, Lx, nblocks, 'single');
ib=0;
for iy = 1:numBlocks(1)
    if numBlocks(2) > 1
        for ix = 1:numBlocks(2)
            ib=ib+1;
            gausy = exp(-([1:Ly]' - yB(iy)).^2 / (2*sT(1).^2));
            gausx = exp(-([1:Lx]' - xB(ix)).^2 / (2*sT(2).^2));
            xyMask(:, :, ib) = gausy * gausx';
        end
    else
        ib=ib+1;
        gausy = exp(-([1:Ly]' - yB(iy)).^2 / (2*sT(1).^2));
        xyMask(:, :, ib) = repmat(gausy, 1, Lx);
    end
end

xyMask = xyMask./repmat(sum(xyMask, 3), 1, 1, size(xyMask, 3));
xyMask = reshape(xyMask, Ly*Lx, nblocks);

ops.xyMask    = xyMask;

fprintf('--- using %d blocks in Y\n', ops.numBlocks(1));
fprintf('--- %d pixels/block; avg pixel overlap = %d pixels\n', ...
    round(ops.blockFrac(1)*Ly), ops.pixoverlap(1) );

if ops.numBlocks(2) > 1
    fprintf('--- using %d blocks in X\n', ops.numBlocks(2));
    fprintf('--- %d pixels/block; avg pixel overlap = %d pixels\n', ...
        round(ops.blockFrac(2)*Lx),  ops.pixoverlap(2));
end
