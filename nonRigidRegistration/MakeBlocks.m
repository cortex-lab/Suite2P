function ops = MakeBlocks(ops0)

ops = ops0;
Lx = ops.Lx;
Ly = ops.Ly;
numBlocks = ops.numBlocks;
bpix      = ops.blockPixels;

yB        = linspace(0, 512, numBlocks+1);
yB        = round((yB(1:end-1) + yB(2:end)) / 2);

for iy = 1:numBlocks
    if iy == numBlocks
        yB(iy)  = Ly - floor(bpix/2);
    elseif iy == 1
        yB(iy)  = floor(bpix/2) + 1;
    end
    ops.yBL{iy} = [max(1,yB(iy)-floor(bpix/2)) : ...
        min(Ly,yB(iy)+floor(bpix/2))];
    ops.xBL{iy} = [1:Lx];
end

xyMask = zeros(Ly, Lx, numBlocks, 'single');
for i = 1:numBlocks
    msk = zeros(Ly, Lx, 'single');
    msk(ops.yBL{i},ops.xBL{i}) = 1;
    sT = floor((numBlocks * bpix - Ly) / (numBlocks - 1));
    sT = max(10, sT);
    msk = my_conv(my_conv(msk, sT)',sT)'; 
    xyMask(:,:,i) = msk;
end

xyMask = xyMask./repmat(sum(xyMask, 3), 1, 1, size(xyMask, 3));
xyMask = reshape(xyMask, Ly*Lx, size(xyMask, 3));

ops.xyMask    = xyMask;