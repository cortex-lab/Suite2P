function ops = MakeBlocks(ops)

Lx = ops.Lx;
Ly = ops.Ly;
numBlocks = ops.numBlocks;
bfrac     = 1/(numBlocks-1);
ops.blockFrac = getOr(ops, {'blockFrac'}, bfrac);
bpix      = round(ops.blockFrac * Ly);
sT = floor((numBlocks * bpix - Ly) / (numBlocks - 1));
ops.smoothBlocks = getOr(ops, {'smoothBlocks'}, sT);
sT = ops.smoothBlocks;

yB        = linspace(0, Ly, numBlocks+1);
yB        = round((yB(1:end-1) + yB(2:end)) / 2);

pixoverlap   = zeros(Ly,1);
for iy = 1:numBlocks
    if iy == numBlocks
        yB(iy)  = Ly - floor(bpix/2);
    elseif iy == 1
        yB(iy)  = floor(bpix/2) + 1;
    end
    ops.yBL{iy} = [max(1,yB(iy)-floor(bpix/2)) : ...
        min(Ly,yB(iy)+floor(bpix/2))];
    ops.xBL{iy} = [1:Lx];
    
    pixoverlap(ops.yBL{iy}) = pixoverlap(ops.yBL{iy}) + 1; 
end

pixoverlap = 2*sum(pixoverlap > 1) / numBlocks;
ops.pixoverlap = pixoverlap;

xyMask = zeros(Ly, Lx, numBlocks, 'single');
for i = 1:numBlocks
    msk = zeros(Ly, Lx, 'single');
    msk(ops.yBL{i},ops.xBL{i}) = 1;
    sT = max(10, sT);
    msk = my_conv(my_conv(msk, sT)',sT)'; 
    xyMask(:,:,i) = msk;
end

xyMask = xyMask./repmat(sum(xyMask, 3), 1, 1, size(xyMask, 3));
xyMask = reshape(xyMask, Ly*Lx, size(xyMask, 3));

ops.xyMask    = xyMask;