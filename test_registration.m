function err = test_registration(frame, nTrans, useGPU)
%test_registration Test registration by registering a randomly shifted frame
%   Shifts a frame with nTrans random translations then registers adds
%   nosie to each frame and attempts to register them back to the frame.
%   Returns the mean square error of the registration offsets with the
%   original translations
if nargin < 2
  nTrans = 1000;
end
if nargin < 3
  useGPU = false;
end
frame = single(frame);
dvrSigma =  2;
noiseSigma = std(frame(:))/3;
dvr = [0 0; dvrSigma*randn(nTrans - 1, 2)];

ops = struct('useGPU', useGPU, 'PhaseCorrelation', true, 'mimg', [],...
  'SubPixel', inf, 'registrationUpsample', 1, 'regPrecision', 'same');

movie = repmat(frame, [1 1 nTrans]);

[transFrame, validIdx] = translate_movie(movie, dvr, ops);
transFrame = transFrame(validIdx{1},validIdx{2},:);

transFrame = transFrame + noiseSigma*randn(size(transFrame), 'single');
ops.mimg = transFrame(:,:,1);
%%
tic
nreps = 1;
for i = 1:nreps
  reg_new = registration_offsets(transFrame, ops, false);
end
dt_new = toc;
tic
for i = 1:nreps
  reg_old = registration_offsets_old(transFrame, ops, false);
end
dt_old = toc;
err_old = mean((dvr(:) - reg_old(:)).^2);
err_new = mean((dvr(:) - reg_new(:)).^2);
oldTimesMoreError = (err_old)/(err_new);
fprintf('dt_new=%.2fs, dt_old=%.2fs old gave %.1f times more MSE\n',dt_new,dt_old,oldTimesMoreError);
err = err_new;
% ff = 1:nTrans;
% plot(ff, dvr(:,1)-crv(:,1),'b',ff, dvr(:,2)-crv(:,2),'r');
end
