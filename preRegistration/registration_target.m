function [idx, frame] = registration_target(frames, useGPU)

%% Parameters
[ly, lx, nFrames] = size(frames);
data = single(frames);

if nargin < 2
  useGPU = false;
end

%% Prepare common arrays
% Taper mask
[ys, xs] = ndgrid(1:ly, 1:lx);
ys = abs(ys - mean(ys(:)));
xs = abs(xs - mean(xs(:)));
mY      = max(ys(:)) - 4;
mX      = max(xs(:)) - 4;
slope   = 1.2; % was 2

maskMul = single(1./(1 + exp((ys - mY)/slope)) ./(1 + exp((xs - mX)/slope)));
maskOffset = mean(data(:,:,1))*(1 - maskMul);
% Smoothing filter in frequency domain
sigma = 0.76; %.76 with mask
hgx = exp(-(((0:lx-1) - fix(lx/2))/sigma).^2);
hgy = exp(-(((0:ly-1) - fix(ly/2))/sigma).^2);
hg = hgy'*hgx;
fhg = real(fftn(ifftshift(single(hg/sum(hg(:))))));

eps0 = single(1e-20);
if useGPU
  data = gpuArray(data);
%   peakCorrMatrix = zeros(nFrames, nFrames, 'single', 'gpuArray');
% else
%   peakCorrMatrix = zeros(nFrames, nFrames, 'single');
end

corrRanges = zeros(nFrames, 3, 'single');

[a, b] = ndgrid(1:nFrames, 1:nFrames);
data = fft2(bsxfun(@plus, maskOffset, bsxfun(@times, maskMul, data)));
data = data./(eps0 + abs(data));
cdata = conj(data);
for fi = 1:nFrames
  % compute correlation map of frame fi with every other frame
  cmap = bsxfun(@times, data(:,:,a(:,fi)).*cdata(:,:,b(:,fi)), fhg);
  cmap = real(ifft2(cmap));
  % find peak of each correlation map
  peaks = gather(max(reshape(cmap, ly*lx, nFrames), [], 1));
  % compute first quartile of peak correlations, ignoring result
  % from correlation of this frame with itself.
  corrRanges(fi,:) = prctile(peaks((1:nFrames) ~= fi), [25 50 75]);
end
% find the frame with the best correlations with other frames
[bestCorr, idx] = max(corrRanges(:,1));
% plot(corrRanges)
frame = frames(:,:,idx);
% peakCorrMatrix = gather(peakCorrMatrix);
% [mxi, idx] = max(prctile(peakCorrMatrix, 25, 1))
% [mnx, mni] = min(prctile(peakCorrMatrix, 25, 2));
% figure, imagesc(peakCorrMatrix);
% figure, plot(1:nFrames, peakCorrMatrix(mxi,:), 1:nFrames, peakCorrMatrix(mni,:))
end