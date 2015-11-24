function [tdata, validIdx] = translate_movie(data, dv, ops)

%% Parameters
[ly, lx, nFrames] = size(data);
if nargin < 3 
  ops = [];
end
subpixel = getOr(ops, {'subPixel' 'SubPixel'}, 1);
if isfinite(subpixel)
  dv = round(subpixel*dv)./subpixel;
end
useGPU = getOr(ops, 'useGPU', false);

dv = permute(dv, [2 3 1]);
fy = ifftshift((-fix(ly/2):ceil(ly/2) - 1)/ly)';% freq along first dimension
fx = ifftshift((-fix(lx/2):ceil(lx/2) - 1)/lx); % freq along second dimension

if useGPU
  batchSize = 32;
  fx = gpuArray(fx);
  fy = gpuArray(fy);
  dv = gpuArray(dv);
else
  batchSize = 3;
end
tdata = zeros([ly, lx, nFrames], 'like', data);
%% Work through data in batches
nBatches = ceil(nFrames/batchSize);
for bi = 1:nBatches
  fi = (bi - 1)*batchSize + 1:min(bi*batchSize, nFrames);
  if useGPU
    batchData = gpuArray(single(data(:,:,fi)));
  else
    batchData = single(data(:,:,fi));
  end
  phaseShift = bsxfun(@times,...
    exp(-1j*2*pi*bsxfun(@times, fy, dv(1,:,fi))),... y rotation
    exp(-1j*2*pi*bsxfun(@times, fx, dv(2,:,fi)))); % x rotation
  tdata(:,:,fi) = gather(real(ifft2(fft2(batchData).*phaseShift)));
end
if nargout > 1
  dvMax = max(0, ceil(max(gather(dv), [], 3)));
  dvMin = min(0, floor(min(gather(dv), [], 3)));
  validX = (1 + dvMax(2)):(lx + dvMin(2));
  validY = (1 + dvMax(1)):(ly + dvMin(1));
  validIdx = {validY validX};
end

end
