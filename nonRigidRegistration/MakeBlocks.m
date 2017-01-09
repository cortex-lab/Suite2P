function ops = MakeBlocks(ops)

Lx = ops.Lx;
Ly = ops.Ly;
numBlocks = ops.numBlocks;
bfrac     = 1/max(2,(numBlocks-3));
ops.blockFrac = getOr(ops, {'blockFrac'}, bfrac);
bpix      = round(ops.blockFrac * Ly);

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

sT        = mean(diff(yB)) * 2/3;
ops.smoothBlocks = getOr(ops, {'smoothBlocks'}, sT);
ops.smoothBlocks = max(10, ops.smoothBlocks);
sT        = ops.smoothBlocks;

xyMask = zeros(Ly, numBlocks, 'single');
for iy = 1:numBlocks
    gaus = exp(-([1:Ly]' - yB(iy)).^2 / (2*sT^2));
    xyMask(:, iy) = gaus;
end
xyMask = xyMask./repmat(sum(xyMask, 2), 1, size(xyMask, 2));

ops.xyMask    = xyMask;