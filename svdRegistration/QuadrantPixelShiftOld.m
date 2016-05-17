function [ops,pixShift] = QuadrantPixelShift(ops,yB,xB,npix,A,B)

ib = 0;
for iy = 1:length(yB)-1
  for ix = 1:length(xB)-1
    ib = ib+1;
    ops.yBL{ib} = (yB(iy)+1):yB(iy+1);
    ops.xBL{ib} = (xB(ix)+1):xB(ix+1);
  end
end
numBlocks = ib;

xyMask = zeros(npix, npix, numBlocks, 'single');
qShift = zeros(numBlocks,2,'single');
for i = 1:numBlocks
    msk = zeros(npix, npix, 'single');
    msk(ops.yBL{i},ops.xBL{i}) = 1;
    sT = 25;
    msk = my_conv(my_conv(msk, sT)',sT)'; 
    xyMask(:,:,i) = msk;
    
    MQ2 = B(ops.yBL{i},ops.xBL{i});
    MQ1 = A(ops.yBL{i},ops.xBL{i});
    [cx,ix] = regZ(MQ1,MQ2);
    qShift(i,:) = ix;
end
xyMask = xyMask./repmat(sum(xyMask, 3), 1, 1, numBlocks);
pixShift = zeros(npix,npix,2,'single');
for i = 1:numBlocks
  for j = 1:2
    pixShift(:,:,j) = pixShift(:,:,j) + xyMask(:,:,i) * qShift(i,j);
  end
end

end