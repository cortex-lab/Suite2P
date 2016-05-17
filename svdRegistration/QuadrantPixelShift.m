function [ops,pixShift,ccMean,ccx,cX] = QuadrantPixelShift(ops,npix,A,B,isGPU,isSVD)

%%%% params TO BE PUT IN OPTIONS
nQuadrants  = 8;
nLocalPix   = 5; % only deviate this far from mean shift!

yB = round(linspace(0,npix,nQuadrants+1));
nWind   = 20;  %%% make sure all quadrants are same size
xB = yB;

if isGPU==1
    A = gpuArray(A);
    B = gpuArray(B);
end

[ops,numBlocks] = MakeQuadrants(ops,yB,xB);


% compute mean shift
if isSVD == 1
    [~,ixAll,ccMean] = SVDPixelShift(A,B,nWind);
else
    [~,ixAll,ccMean] = regZlocal(A,B,nWind);
end

%disp(ixAll);
if isGPU==1
    ccMean = gather(ccMean);
    ixAll = gather(ixAll);
end

xyMask = zeros(npix, npix, numBlocks, 'single');
qShift = zeros(numBlocks,2,'single');
for i = 1:numBlocks
    % compute masks for each quadrant for smoothing
    xyMask(:,:,i) = QuadrantMask(ops,npix,i);

    % compute correlation matrix per quadrant
    if isSVD == 1
       MQ2 = B(ops.yBL{i},ops.xBL{i},:);
       MQ1 = A(ops.yBL{i},ops.xBL{i},:);
       [~,~,cc] = SVDPixelShift(MQ1,MQ2,nWind);
    else
        MQ2 = B(ops.yBL{i},ops.xBL{i});
        MQ1 = A(ops.yBL{i},ops.xBL{i});
        [~,~,cc] = regZlocal(MQ1,MQ2,nWind);
    end
    [ii,jj] = ind2sub([nQuadrants nQuadrants],i);
    ccx(:,:,ii,jj) = cc;
end

% smooth correlation matrix over quadrants
ccx = my_conv2(ccx,2,[3 4]);


if isGPU == 1
    ccx = gather(ccx);
end
for i = 1:numBlocks
    ccx0 = ccx(:,:,i);
    %%%% find best subpixel alignment for each quadrant
    % in nLocalPix region around full image shift
    [cx,ix] = SubPixel2D(ccx0,ixAll,nLocalPix);
    ix = FindRegInds(ix,size(ccx0,1),size(ccx0,2));
    cX(i) = cx;                           
    qShift(i,:) = ix;
end


xyMask = xyMask./repmat(sum(xyMask, 3), 1, 1, numBlocks);
pixShift = zeros(npix,npix,2,'single');
for i = 1:numBlocks
  for j = 1:2
    pixShift(:,:,j) = pixShift(:,:,j) + xyMask(:,:,i) * qShift(i,j);
  end
end