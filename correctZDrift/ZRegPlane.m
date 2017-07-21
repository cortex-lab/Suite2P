
function [cMAX, ix, cZ] = ZRegPlane(m1,m2,indsZ,useGPU)


cMAX = 0; ixMAX = [0 0]; iZ = 0;
if useGPU == 1
    cc0 = gpuArray.zeros(size(m2,1),size(m2,2),length(indsZ),'single');
    cZ  = gpuArray.zeros(length(indsZ),1,'single');
else
    cc0 = zeros(size(m2,1),size(m2,2),length(indsZ),'single');
    cZ  = zeros(length(indsZ),1,'single');
end
k=0;

% compute phase-correlation with each plane in z-stack
for z = indsZ(:)'
    k=k+1;
  cc = real(ifft(ifft(m1(:,:,z) .* conj(m2),[],1),[],2));
  cc0(:,:,k) = cc;
  
  [cmax,ixmax] = max(cc(:));
  [ix1,ix2] = ind2sub(size(cc),ixmax);
  ixmax = [ix1 ix2];
  if cmax > cMAX
      cMAX  = cmax;
      ixMAX = ixmax;
      iZ    = z;
  end
  cZ(k) = cmax;
end

cc0 = fftshift(cc0);
[cmax0,ixmax0] = max(cc0(:));
[ix1,ix2,ix3] = ind2sub(size(cc0),ixmax0);
ix0 = [ix1 ix2 ix3];

%%% sub-pixel alignment in z
% first average in y,x around peak
dZ=[-4:4];
zinds = (unique(min(size(cc0,3),max(1,ix0(3)+dZ))));
zinds(zinds>size(cc0,3) | zinds<1) = [];
ccZ=cc0(ix0(1)+[-1:1],ix0(2)+[-1:1],zinds);
ccZ = squeeze(sum(sum(ccZ,1),2));
zinds = zinds - ix0(3);
if useGPU==1
    ccZ = gather(ccZ);
    zinds=gather(zinds);
    
    cMAX = gather(cMAX);
    ixMAX = gather(ixMAX);
    cZ =gather(cZ);
end
% then interp in z and find max
dZ2 = linspace(-2,2,100);
cZi=interp1(zinds,ccZ,dZ2,'spline');
[cmax,izmax] = max(cZi);
dZi = dZ2(izmax);
iZ = iZ+dZi;

[nX, nY] = size(m2);
ixMAX = FindRegInds(ixMAX,[nX nY]);

ix = [ixMAX iZ];

if useGPU==1
    ix = gather(ix);
end
    