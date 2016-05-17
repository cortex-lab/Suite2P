
if ~exist('loaded');
    clear dat;
    fs = dir('~/DATA/zStacks/SVDs/*.mat');
    for i = [1:6]
        dat{i} = load(sprintf('~/DATA/zStacks/SVDs/%s',fs(i).name));
    end

    loaded = 1;
end

%%
nBest = 4;
clear U;
npix = 512;
Upad = zeros(npix,npix,size(dat{1}.U,3),'single');
Upad(dat{nBest}.ops.yrange,dat{nBest}.ops.xrange,:) = dat{nBest}.U;
dat{nBest}.ops.yrange2=dat{nBest}.ops.yrange;
dat{nBest}.ops.xrange2=dat{nBest}.ops.xrange;
U{nBest} = Upad;
dinds = true(1,length(dat));
dinds(nBest) = 0;
dinds = find(dinds);
for nD = dinds
  % reference image
  A = single(dat{nBest}.ops.mimg1);
  % new image
  B = single(dat{nD}.ops.mimg1);

  % compute quadrant based shifts
  % yB and xB are beginning and ends of quadrants
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
  drawnow;
  % aligned new U to reference image
  %  isGPU = 1;
  U0 = ShiftSVD(dat{nD},pixShift,npix,isGPU);
  U{nD} = U0;

  yrange = dat{nD}.ops.yrange;
  xrange = dat{nD}.ops.xrange;
  ymax = round(max(0,max(max(pixShift(:,:,1)))));
  xmax = round(max(0,max(max(pixShift(:,:,2)))));
  ymin = round(min(0,min(min(pixShift(:,:,1)))));
  xmin = round(min(0,min(min(pixShift(:,:,2)))));
  dat{nD}.ops.yrange2 = [yrange(1)-ymin+1 : yrange(end)-ymax];
  dat{nD}.ops.xrange2 = [xrange(1)-xmin+1 : xrange(end)-xmax];

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

for nD = 1:length(U)
  for i = 1:size(U{nD},3)
    U{nD}(:,:,i) = U{nD}(:,:,i) * dat{nD}.Sv(i)^.5;
  end
end

Uall = U{1}(yrange,xrange,:);
for nD = 2:length(dat)
    Uall = cat(3,Uall,U{nD}(yrange,xrange,:));
end

%% take top 1000 components
Uall = reshape(Uall,[],size(Uall,ndims(Uall)));
[Un Sv Vn] = svdecon(Uall'*Uall);

%%
VF = normc(Uall * Un);
UDset = U;
Sv = diag(Sv)/(nD/2);

[ops,stat,res] = fast_clustering_with_neuropil(ops,VF,Sv);


B = single(dat{nBest}.ops.mimg1);
for nD = dinds
    A = single(dat{nD}.ops.mimg1);
    [res0,pixInv] = ShiftMasks(res,ops,A,B,isGPU)
    dat{nD}.res    = res0;
    dat{nD}.pixInv = pixInv;
end












