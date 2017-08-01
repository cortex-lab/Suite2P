% needs ops.RegFile, ops.xrange, ops.yrange, ops.NavgFramesSVD
% Ly, Lx are the size of each frame. nimgbatch: the number of frames loaded per batch
% nt0 is the number of timepoints to bin over. If a sixth argument is
% present, it does not subtract the mean of each batch. 
function mov = loadAndBin(ops, Ly, Lx, nimgbatch, nt0, clustModel)

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
    badi = ops.badframes(ix + [1:size(data,3)]);
    data(:,:, badi) = [];
    
    % subtract off the mean of this batch
    if nargin<=5
        data = bsxfun(@minus, data, mean(data,3));
    end
    %     data = bsxfun(@minus, data, ops.mimg1);
    
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

mov = mov(:, :, 1:ix);
