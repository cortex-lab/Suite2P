% compute SVD of data and save to file
function [ops, U, Sv] = get_svdcomps(ops)

% load(sprintf('%s/%s/%s/regops_%s_%s_plane%d.mat', ops.ResultsSavePath, ops.mouse_name, ops.date, ...
%     ops.mouse_name, ops.date, ops.iplane))

iplane = ops.iplane;

[Ly, Lx] = size(ops.mimg);

ntotframes          = ceil(sum(ops.Nframes));
% number of frames used to compute SVD
ops.NavgFramesSVD   = min(ops.NavgFramesSVD, ntotframes);
% size of binning (in time)
nt0 = ceil(ntotframes / ops.NavgFramesSVD);

if isfield(ops, 'chunk_align') && ~isempty(ops.chunk_align); chunk_align   = ops.chunk_align(iplane);
else chunk_align = 1; end

if chunk_align>9
    nt0 =  ops.chunk_align;
end
ops.NavgFramesSVD = floor(ntotframes/nt0);
nimgbatch = nt0 * floor(2000/nt0);

%% load the data 

ix = 0;
fid = fopen(ops.RegFile, 'r');
mov = zeros(numel(ops.yrange), numel(ops.xrange), ops.NavgFramesSVD, 'single');

while 1
    data = fread(fid,  Ly*Lx*nimgbatch, '*int16');
    if isempty(data)
        break;
    end
    data = reshape(data, Ly, Lx, []);
    
    % subtract off the mean of this batch
%     data = data - repmat(ops.mimg1, 1, 1, size(data,3));
    
    nSlices = nt0*floor(size(data,3)/nt0);
    if nSlices ~= size(data,3)
        data = data(:,:, 1:nSlices);
    end
    
    % bin data
    data = reshape(data, Ly, Lx, nt0, []);
    data = single(data);
    davg = squeeze(mean(data,3));
    
    mov(:,:,ix + (1:size(davg,3))) = davg(ops.yrange, ops.xrange, :);
    
    ix = ix + size(davg,3);
end
fclose(fid);
mov = mov(:, :, 1:ix);

%% SVD options

% number of SVD components kept
ops.nSVD = min(ops.nSVD, size(mov,3));
%
mov             = reshape(mov, [], size(mov,3));
% mov             = mov./repmat(mean(mov.^2,2).^.5, 1, size(mov,2));

% compute covariance matrix of frames
if ops.useGPU
    COV             = gpuBlockXtX(mov)/size(mov, 1);
else
    COV             = mov' * mov/size(mov,1);
end

ops.nSVD = min(size(COV,1)-2, ops.nSVD);

% take SVD of covariance matrix and keep ops.nSVD components
if ops.nSVD<1000 || size(COV,1)>1e4
    [V, Sv]          = eigs(double(COV), ops.nSVD);
else
    if ops.useGPU
        gpuCOV = gpuArray(double(COV));
        [V, Sv]         = svd(gpuCOV);
        clear gpuCOV;
        V = gather(single(V));
        Sv = gather(single(Sv));
    else
         [V, Sv]         = svd(COV);
    end
    V               = V(:, 1:ops.nSVD);
    Sv              = Sv(1:ops.nSVD, 1:ops.nSVD);
end
%%

% compute U (normalized spatial masks... pixels x components)
if ops.useGPU
    U               = normc(gpuBlockXY(mov, V));
else
    U               = normc(mov*V);
end
U               = single(U);
Sv              = single(diag(Sv));

if ~exist(ops.ResultsSavePath, 'dir')
    mkdir(ops.ResultsSavePath)
end

% project spatial masks onto raw data
fid = fopen(ops.RegFile, 'r');
ix = 0;
Fs = zeros(ops.nSVD, sum(ops.Nframes), 'single');
while 1
    data = fread(fid,  Ly*Lx*nimgbatch, '*int16');
    if isempty(data)
        break;
    end
    data = reshape(data, Ly, Lx, []);
    
    % subtract off the mean of this batch
    %         data = data - repmat(mean(data,3), 1, 1, size(data,3));
%     data = data - repmat(ops.mimg1, 1, 1, size(data,3));
    data = data(ops.yrange, ops.xrange, :);
    data = single(data);
    if ops.useGPU
        Fs(:, ix + (1:size(data,3))) = gpuBlockXtY(U, reshape(data, [], size(data,3)));
    else
        Fs(:, ix + (1:size(data,3))) = U' * reshape(data, [], size(data,3));
    end
    
    ix = ix + size(data,3);
end
fclose(fid);

if ~exist(ops.ResultsSavePath, 'dir')
    mkdir(ops.ResultsSavePath);
end

totF = [0 cumsum(ops.Nframes)];
for iexp = 1:length(ops.expts)
    Vcell{iexp} = Fs(:, (1+ totF(iexp)):totF(iexp+1));
end

%%% save SVDs
U = reshape(U, numel(ops.yrange), numel(ops.xrange), []);
try % this is faster, but is limited to 2GB files
    save(sprintf('%s/SVD_%s_%s_plane%d.mat', ops.ResultsSavePath, ...
        ops.mouse_name, ops.date, iplane), 'U', 'Sv', 'Vcell', 'ops', '-v6');
catch % this takes a bit less space, but is significantly slower
    save(sprintf('%s/SVD_%s_%s_plane%d.mat', ops.ResultsSavePath, ...
        ops.mouse_name, ops.date, iplane), 'U', 'Sv', 'Vcell', 'ops');
end

% keyboard;