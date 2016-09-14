function msk = QuadrantMask(yBL,xBL,nY,nX,i)

msk = zeros(nY, nX, 'single');
msk(yBL{i},xBL{i}) = 1;
sT = 25;
msk = my_conv2(msk,sT,[1 2]);