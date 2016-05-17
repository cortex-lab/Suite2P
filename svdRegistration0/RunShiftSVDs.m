function [Uall,dat,ops] = RunShiftSVDs(Uall,dat,nBest,isGPU)

npix = size(Uall,1);

yrange0 = dat{nBest}.ops.yrange;
xrange0 = dat{nBest}.ops.xrange;
dat{nBest}.ops.yrange2=yrange0;
dat{nBest}.ops.xrange2=xrange0;
dinds = true(1,length(dat));
dinds(nBest) = 0;
dinds = find(dinds);
% reference image
A = single(dat{nBest}.ops.mimg1);
for nD = dinds
  % new image
  B = single(dat{nD}.ops.mimg1);

  % compute quadrant based shifts
  % yB and xB are beginning and ends of quadrants
  % only allow 5 x 5 quadrants
  yB = cumsum([0 171 170 171]);  %[0 171 170 171]
  xB = cumsum([0 171 170 171]);

  [ops,pixShift] = QuadrantPixelShift(dat{nD}.ops,yB,xB,npix,A,B);
  disp(max(max(pixShift(:,:,1))))
  disp(max(max(pixShift(:,:,2))))

  dat{nD}.ops = ops;
  dat{nD}.pixShift = pixShift;

  % check that quadrant shifting helped
  % shift B to A
  B0 = ShiftMeanImage(A,B,pixShift,npix);
%   drawnow;
  % aligned new U to reference image
  %  isGPU = 1;
  U0 = ShiftSVD(Uall(:,:,:,nD), pixShift,npix,isGPU);
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

dat{nBest}.ops.yrange2 = yrange;
dat{nBest}.ops.xrange2 = xrange;
for nD = dinds
	dat{nD}.ops.yrange2 = yrange0(1) - dat{nD}.ops.yrange(1) + yrange;
	dat{nD}.ops.xrange2 = xrange0(1) - dat{nD}.ops.xrange(1) + xrange;
end

for nD = 1:size(Uall,4)
  for i = 1:size(Uall,3)
    Uall(:,:,i,nD) = Uall(:,:,i,nD) * dat{nD}.Sv(i)^.5;
  end
end

% Uall = Uall(ops.yrange,ops.xrange,:,:);