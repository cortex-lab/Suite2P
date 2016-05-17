

function [res0,pixShift] = ShiftMasks(res,ops,A,B,isGPU)

npix = 512;
%nY = length(ops.yrange);
%nX = length(ops.xrange);
nY = res.Ly;
nX = res.Lx;
res0 = res;

yB = cumsum([0 171 170 171]);  
xB = cumsum([0 171 170 171]);

[~,pixShift] = QuadrantPixelShift(ops,yB,xB,npix,A,B);


for k=1:3
    switch k
      case 1
        U1 = res.M;
        U1 = reshape(U1,nY,nX);
      case 2
        U1 = res.S;
        U1 = reshape(U1,[nY,nX,size(U1,2)]);
      case 3
        U1 = res.iclust;
        U1 = reshape(U1,nY,nX);
    end

    nComps = size(U1,3);
    U0 = zeros(npix,npix,nComps,'single');
    U0(ops.yrange,ops.xrange,:) = U1;

    x=[1:npix]; y=[1:npix];
    [xx0,yy0] = meshgrid(x,y);
    xx = xx0 - pixShift(:,:,2);
    yy = yy0 - pixShift(:,:,1);

    if isGPU == 1
        x = gpuArray(single(x));
        y = gpuArray(single(y));
        xx = gpuArray(single(xx));
        yy = gpuArray(single(yy));
    end

    for i = 1:nComps
        if isGPU == 1
            Un = gpuArray(U0(:,:,i));
        else
            Un = U0(:,:,i);
        end
        Un = interp2(x,y,Un,xx,yy,'nearest');
        % nan's tell new xrange and yrange
        if k==1
            yv = gather(Un(:,250));
            xv = gather(Un(250,:));
            res0.yrange = find(yv~=0 & ~isnan(yv));
            res0.xrange = find(xv~=0 & ~isnan(xv));
            res0.xrange=res0.xrange(:);
        end
        Un(isnan(Un(:))) = 0;
        if isGPU == 1
            U0(:,:,i) = gather(Un);
        else
            U0(:,:,i) = Un;
        end
    end

    switch k
      case 1
        res0.M = U0(res0.yrange, res0.xrange);
        res0.M = res0.M(:)';
      case 2
        res0.S = reshape(U0(res0.yrange, res0.xrange, :), [],size(U0,3));
      case 3
        res0.iclust = U0(res0.yrange, res0.xrange);
        res0.iclust = res0.iclust(:)';
    end
end

res0.Ly = length(res0.yrange);
res0.Lx = length(res0.xrange);