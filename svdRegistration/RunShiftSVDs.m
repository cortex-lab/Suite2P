function [Uall,dat,ops] = RunShiftSVDs(Uall,dat,nBest,isGPU,isSVD)

npix = size(Uall,1);

yrange0 = dat{nBest}.ops.yrange;
xrange0 = dat{nBest}.ops.xrange;
dat{nBest}.ops.yrange2=yrange0;
dat{nBest}.ops.xrange2=xrange0;
dinds = true(1,length(dat));
dinds(nBest) = 0;
dinds = find(dinds);
% reference image
if isSVD == 1
    A = single(Uall(:,:,1:10,nBest));
else
    A = single(dat{nBest}.ops.mimg1);
end

for nD = dinds
  % new image
  if isSVD == 1
      B = single(Uall(:,:,1:10,nD));
  else
      B = single(dat{nD}.ops.mimg1);
  end

  % compute quadrant based shifts
  [ops,pixShift,ccMean,ccQuad,cX] = QuadrantPixelShift(dat{nD}.ops, ...
                                                    npix,A,B,isGPU,isSVD);

  %clf;
  %for j = 1:2
  %    subplot(1,2,j),
  %    imagesc(pixShift(:,:,j));
  %    colorbar;
  %end
  %drawnow;
  
  
  %disp(max(max(abs(pixShift(:,:,1)))))
  %disp(max(max(abs(pixShift(:,:,2)))))


  dat{nD}.ops = ops;
  dat{nD}.pixShift = pixShift;

  % check that quadrant shifting helped
  % shift B to A
  % aligned new U to reference image
  %  isGPU = 1;
  U0 = ShiftSVD(Uall(:,:,:,nD), pixShift,npix,isGPU);
  %[U1,cX] = SVDProj(A,U0(:,:,1:10),isGPU);

  Uall(:,:,:,nD) = U0;
  clear U0;

  yrange = dat{nD}.ops.yrange;
  xrange = dat{nD}.ops.xrange;
  ymax = round(max(0,max(max(pixShift(:,:,1)))));
  xmax = round(max(0,max(max(pixShift(:,:,2)))));
  ymin = round(min(0,min(min(pixShift(:,:,1)))));
  xmin = round(min(0,min(min(pixShift(:,:,2)))));
  dat{nD}.ops.yrange2 = [yrange(1)+ymax : yrange(end)+ymin];
  dat{nD}.ops.xrange2 = [xrange(1)+xmax : xrange(end)+xmin];

  %if nD == 2
  %    keyboard;
  %end
end

%% combine all U's
yrange = dat{1}.ops.yrange2;
xrange = dat{1}.ops.xrange2;

for nD = 2:length(dat)
    yrange = intersect(yrange,dat{nD}.ops.yrange2);
    xrange = intersect(xrange,dat{nD}.ops.xrange2);
end
ops.yrange = yrange;
ops.xrange = xrange;

dat{nBest}.ops.yrange0 = yrange;
dat{nBest}.ops.xrange0 = xrange;
for nD = dinds
    dat{nD}.ops.yrange0 = yrange0(1) - dat{nD}.ops.yrange(1) + yrange;
    dat{nD}.ops.xrange0 = xrange0(1) - dat{nD}.ops.xrange(1) + xrange;
end

for nD = 1:size(Uall,4)
  for i = 1:size(Uall,3)
    Uall(:,:,i,nD) = Uall(:,:,i,nD) * dat{nD}.Sv(i)^.5;
  end
end
