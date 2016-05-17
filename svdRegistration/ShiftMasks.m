function [res0,pixShift] = ShiftMasks(res,opsSingle,opsMulti,A,B, ...
                                      isGPU,isSVD)

ops0 = opsMulti;
ops  = opsSingle;

npix = 512;
nY = res.Ly;
nX = res.Lx;
res0 = res;

[~,pixShift,cX] = QuadrantPixelShift(ops,npix,A,B,isGPU,isSVD);
fprintf('yshift %2.2f xshift %2.2f\n',...
	mean(mean(pixShift(:,:,1))),mean(mean(pixShift(:,:,2))));
for k=1:3
    switch k
      case 1
        U1 = res.M;
        U1 = reshape(U1,nY,nX,[]);
      case 2
        U1 = res.S;
        U1 = reshape(U1,[nY,nX,size(U1,2)]);
      case 3
        U1 = res.iclust;
        U1 = reshape(U1,nY,nX);
    end

    nComps = size(U1,3);
    U0 = zeros(npix,npix,nComps,'single');
    U0(ops0.yrange,ops0.xrange,:) = U1;
    clear U1;

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
            xrange  = ops.xrange;
            yrange  = ops.yrange;
        end
        Un(isnan(Un(:))) = 0;
        if isGPU == 1
            U0(:,:,i) = gather(Un);
        else
            U0(:,:,i) = Un;
        end
    end

    %check pixel shifts
    %[cX,ix] = CheckPixelShift(A,U0,isGPU);

    switch k
      case 1
        res0.M = U0(yrange,xrange);
        res0.M = res0.M(:)';
      case 2
        res0.S = reshape(U0(yrange, xrange, :), [],size(U0,3));
        res0.S = normc(res0.S);
      case 3
        res0.iclust = U0(yrange, xrange);
        res0.iclust = res0.iclust(:);
    end
end

res0.Ly = length(yrange);
res0.Lx = length(xrange);

