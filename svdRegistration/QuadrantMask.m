function msk = QuadrantMask(ops,npix,i)

msk = zeros(npix, npix, 'single');
msk(ops.yBL{i},ops.xBL{i}) = 1;
sT = 25;
msk = my_conv2(msk,sT,[1 2]);