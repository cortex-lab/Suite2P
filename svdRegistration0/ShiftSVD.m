function U = ShiftSVD(U,pixShift,npix,isGPU)

nComps = size(U,3);
npix = size(U,1);
% don't declare on GPU -- may be too big

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
    Un = gpuArray(U(:,:,i));
  else
    Un = U(:,:,i);
  end
  Un = interp2(x,y,Un,xx,yy);
  Un(isnan(Un(:))) = 0;
  if isGPU == 1
    U(:,:,i) = gather(Un);
  else
    U(:,:,i) = Un;
  end
end

end