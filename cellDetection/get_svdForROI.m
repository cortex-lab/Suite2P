% compute SVD and project onto normalized data
% save if ops.writeSVDroi
function [ops, U, U2] = get_svdForROI(ops)


sig = ops.diameter/10.;
[Ly, Lx] = size(ops.mimg1);

ntotframes          = ceil(sum(ops.Nframes));
% number of frames to use for SVD
bin_min = round(ntotframes/ops.NavgFramesSVD);
bin_min = max(1, bin_min);
bin_tau = round(ops.tau * ops.fs);
nt0 = max(bin_min, bin_tau);
ops.NavgFramesSVD = floor(ntotframes/nt0);

nimgbatch = nt0 * 100;
nimgbatch = min(ntotframes, nimgbatch);
nimgbatch = nt0 * floor(nimgbatch/nt0);
disp(nt0)

% ops.NavgFramesSVD   = min(ops.NavgFramesSVD, ntotframes);
% nt0 = ceil(ntotframes / ops.NavgFramesSVD);
% ops.NavgFramesSVD = floor(ntotframes/nt0);
% nimgbatch = nt0 * floor(2000/nt0);
%%
mov = loadAndBin(ops, Ly, Lx, nimgbatch, nt0);

% SVD options
ops.nSVDforROI = min(ops.nSVDforROI, round(size(mov,3)/2));

% smooth spatially to get high SNR SVD components
for i = 1:size(mov,3)
    I = mov(:,:,i);
    I = my_conv2(I, sig, [1 2]); %my_conv(my_conv(I',ops.sig)', ops.sig);
    mov(:,:,i) = I;
end
% for i = 1:size(mov,2)
%     mov(:,i,:) = mov(:,i,:) - my_conv2(mov(:,i,:), 1,3);
% end

mov             = reshape(mov, [], size(mov,3));
% compute noise variance across frames (assumes slow signal)
sdmov           = mean(diff(mov, 1, 2).^2, 2).^.5;
sdmov           = reshape(sdmov, numel(ops.yrange), numel(ops.xrange));
sdmov           = max(1e-10, sdmov);
ops.sdmov       = sdmov;

% normalize pixels by noise variance
mov             = bsxfun(@rdivide, mov, sdmov(:));

% compute covariance of frames
COV             = mov' * mov/size(mov,1);

% compute SVD of covariance matrix
[V, Sv]          = eigs(double(COV), ops.nSVDforROI);
% Sv              = single(diag(Sv));

U               = mov * V;
U               = single(U);
% reshape U to frame size
U = reshape(U, numel(ops.yrange), numel(ops.xrange), []);

% compute spatial masks (U = mov * V)
mov = loadAndBin(ops, Ly, Lx, nimgbatch, nt0);
if getOr(ops, 'smooth_masks', 1)
    % smooth spatially to get high SNR SVD components
    for i = 1:size(mov,3)
        I = mov(:,:,i);
        I = my_conv2(I, sig/2, [1 2]); %my_conv(my_conv(I',ops.sig)', ops.sig);
        mov(:,:,i) = I;
    end
    for i = 1:size(mov,2)
%         mov(:,i,:) = mov(:,i,:) - my_conv2(mov(:,i,:), 1,3);
    end
end
mov             = reshape(mov, [], size(mov,3));
% compute noise variance across frames (assumes slow signal)
sdmov           = mean(diff(mov, 1, 2).^2, 2).^.5;
sdmov           = reshape(sdmov, numel(ops.yrange), numel(ops.xrange));
sdmov           = max(1e-10, sdmov);
ops.sdmov       = sdmov;

% normalize pixels by noise variance
mov             = bsxfun(@rdivide, mov, sdmov(:));

U2               = mov * V;
U2               = single(U2);
% reshape U to frame size
U2 = reshape(U2, numel(ops.yrange), numel(ops.xrange), []);

end

% needs ops.RegFile, ops.xrange, ops.yrange, ops.NavgFramesSVD
% Ly, Lx are the size of each frame. nimgbatch: the number of frames loaded per batch
% nt0 is the number of timepoints to bin over. If a sixth argument is
% present, it does not subtract the mean of each batch. 
function mov = loadAndBin(ops, Ly, Lx, nimgbatch, nt0)

ix = 0;
fid = fopen(ops.RegFile, 'r');
mov = zeros(numel(ops.yrange), numel(ops.xrange), ops.NavgFramesSVD, 'single');
ij = 0;
while 1
    % load frames
    data = fread(fid,  Ly*Lx*nimgbatch, '*int16');
    if isempty(data)
        break;
    end
    data = single(data);
    data = reshape(data, Ly, Lx, []);
    
    % ignore bad frames
%     badi = ops.badframes(ix + [1:size(data,3)]);
%     data(:,:, badi) = [];
    
    % subtract off the mean of this batch    
    data = bsxfun(@minus, data, mean(data,3));
    
    nSlices = nt0*floor(size(data,3)/nt0);
    if nSlices~=size(data,3)
        data = data(:,:, 1:nSlices);
    end
    
    % bin data
    data = reshape(data, Ly, Lx, nt0, []);
    davg = squeeze(mean(data,3));
    
    mov(:,:,ix + (1:size(davg,3))) = davg(ops.yrange, ops.xrange, :);
    
    ix = ix + size(davg,3);
    ij = ij + 1;
end
fclose(fid);
end

