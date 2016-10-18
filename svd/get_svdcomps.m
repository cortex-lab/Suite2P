function [ops, U, Sv] = get_svdcomps(ops)

% load(sprintf('%s/%s/%s/regops_%s_%s_plane%d.mat', ops.ResultsSavePath, ops.mouse_name, ops.date, ...
%     ops.mouse_name, ops.date, ops.iplane))

iplane = ops.iplane;

[Ly, Lx] = size(ops.mimg);

ntotframes          = ceil(sum(ops.Nframes));
ops.NavgFramesSVD   = min(ops.NavgFramesSVD, ntotframes);
nt0 = ceil(ntotframes / ops.NavgFramesSVD);

if isfield(ops, 'chunk_align') && ~isempty(ops.chunk_align); chunk_align   = ops.chunk_align(iplane);
else chunk_align = 1; end

if chunk_align>9
    nt0 =  ops.chunk_align;
end
ops.NavgFramesSVD = floor(ntotframes/nt0);
nimgbatch = nt0 * floor(2000/nt0);

ix = 0;
fid = fopen(ops.RegFile, 'r');

mov = zeros(Ly, Lx, ops.NavgFramesSVD, 'single');

while 1
    data = fread(fid,  Ly*Lx*nimgbatch, '*int16');
    if isempty(data)
        break;
    end
    data = single(data);
    data = reshape(data, Ly, Lx, []);
    
    % subtract off the mean of this batch
%     data = data - repmat(ops.mimg1, 1, 1, size(data,3));
    
    irange = 1:nt0*floor(size(data,3)/nt0);
    data = data(:,:, irange);
    
    data = reshape(data, Ly, Lx, nt0, []);
    davg = single(squeeze(mean(data,3)));
    
    mov(:,:,ix + (1:size(davg,3))) = davg;
    
    ix = ix + size(davg,3);
end
fclose(fid);
%%
mov(:, :, (ix+1):end) = [];

mov = mov(ops.yrange, ops.xrange, :);
%% SVD options

ops.nSVD = min(ops.nSVD, size(mov,3));
%
mov             = reshape(mov, [], size(mov,3));
% mov             = mov./repmat(mean(mov.^2,2).^.5, 1, size(mov,2));
COV             = mov' * mov/size(mov,1);

ops.nSVD = min(size(COV,1)-2, ops.nSVD);



if ops.nSVD<1000 || size(COV,1)>1e4
    [V, Sv]          = eigs(double(COV), ops.nSVD);
else
    if ops.useGPU
        [V, Sv]         = svd(gpuArray(double(COV)));
        V = gather(single(V));
        Sv = gather(single(Sv));
    else
         [V, Sv]         = svd(COV);
    end
    V               = V(:, 1:ops.nSVD);
    Sv              = Sv(1:ops.nSVD, 1:ops.nSVD);
end
%%
U               = normc(mov * V);
U               = single(U);
Sv              = single(diag(Sv));

if ~exist(ops.ResultsSavePath, 'dir')
    mkdir(ops.ResultsSavePath)
end

fid = fopen(ops.RegFile, 'r');

ix = 0;
Fs = zeros(ops.nSVD, sum(ops.Nframes), 'single');
while 1
    data = fread(fid,  Ly*Lx*nimgbatch, '*int16');
    if isempty(data)
        break;
    end
    data = single(data);
    data = reshape(data, Ly, Lx, []);
    
    % subtract off the mean of this batch
    %         data = data - repmat(mean(data,3), 1, 1, size(data,3));
%     data = data - repmat(ops.mimg1, 1, 1, size(data,3));
    data = data(ops.yrange, ops.xrange, :);
    Fs(:, ix + (1:size(data,3))) = U' * reshape(data, [], size(data,3));
    
    ix = ix + size(data,3);
end
fclose(fid);

if ~exist(ops.ResultsSavePath, 'dir')
    mkdir(ops.ResultsSavePath);
end

totF = [0 cumsum(ops.Nframes)];
for iexp = 1:length(ops.expts)
    Vcell{iexp} = Fs(:, (1+ totF(iexp)):totF(iexp+1));
    %         save(sprintf('%s/%s/%s/SVD_%s_%s_exp%d_plane%d.mat', ops.ResultsSavePath, ops.mouse_name, ops.date,...
    %             ops.mouse_name, ops.date, ops.expts(iexp), iplane), 'F')
end

U = reshape(U, numel(ops.yrange), numel(ops.xrange), []);
save(sprintf('%s/SVD_%s_%s_plane%d.mat', ops.ResultsSavePath, ...
    ops.mouse_name, ops.date, iplane), 'U', 'Sv', 'Vcell', 'ops');
% keyboard;