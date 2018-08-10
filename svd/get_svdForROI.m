% compute SVD and project onto normalized data
% save if ops.writeSVDroi
function [ops, U, model, U2] = get_svdForROI(ops, clustModel, sm)

% iplane = ops.iplane;
U = []; Sv = []; V = []; Fs = []; sdmov = [];

[Ly, Lx] = size(ops.mimg1);

ntotframes          = ceil(sum(ops.Nframes));
% number of frames to use for SVD
ops.NavgFramesSVD   = min(ops.NavgFramesSVD, ntotframes);
% size of binning (in time)
nt0 = ceil(ntotframes / ops.NavgFramesSVD);


ops.NavgFramesSVD = floor(ntotframes/nt0);
nimgbatch = nt0 * floor(2000/nt0);

mov = loadAndBin(ops, Ly, Lx, nimgbatch, nt0);

% mov = mov - repmat(mean(mov,3), 1, 1, size(mov,3));
%% SVD options
if nargin==1 || ~strcmp(clustModel, 'CNMF')
    ops.nSVDforROI = min(ops.nSVDforROI, size(mov,3));

    % smooth spatially to get high SNR SVD components
    if ops.sig>0.05
        for i = 1:size(mov,3)
            I = mov(:,:,i);
            I = my_conv2(I, ops.sig, [1 2]); %my_conv(my_conv(I',ops.sig)', ops.sig);
            mov(:,:,i) = I;
        end
    end

    mov             = reshape(mov, [], size(mov,3));
    % compute noise variance across frames (assumes slow signal)
    if 1
        sdmov           = mean(diff(mov, 1, 2).^2, 2).^.5;
    else
        sdmov           = mean(mov.^2,2).^.5;
    end
    sdmov           = reshape(sdmov, numel(ops.yrange), numel(ops.xrange));
    sdmov           = max(1e-10, sdmov);
    ops.sdmov       = sdmov;

    % normalize pixels by noise variance
    mov             = bsxfun(@rdivide, mov, sdmov(:));
    model.sdmov     = sdmov;

    % compute covariance of frames
    if ops.useGPU
        COV             = gpuBlockXtX(mov)/size(mov,1);
    else
        COV             = mov' * mov/size(mov,1);
    end
    ops.nSVDforROI = min(size(COV,1)-2, ops.nSVDforROI);

    % compute SVD of covariance matrix
    if ops.useGPU && size(COV,1)<1.2e4
        [V, Sv, ~]      = svd(gpuArray(double(COV)));
        V               = single(V(:, 1:ops.nSVDforROI));
        Sv              = single(diag(Sv));
        Sv              = Sv(1:ops.nSVDforROI);
        %
        Sv = gather(Sv);
        V = gather(V);
    else
        [V, Sv]          = eigs(double(COV), ops.nSVDforROI);
        Sv              = single(diag(Sv));
    end


    if ops.useGPU
        U               = gpuBlockXY(mov, V);
    else
        U               = mov * V;
    end
    U               = single(U);
    % reshape U to frame size
    U = reshape(U, numel(ops.yrange), numel(ops.xrange), []);


    % compute spatial masks (U = mov * V)
    mov = loadAndBin(ops, Ly, Lx, nimgbatch, nt0);
    mov             = reshape(mov, [], size(mov,3));

    if ops.useGPU
        U2               = gpuBlockXY(mov, V);
    else
        U2               = mov * V;
    end
    U2               = single(U2);
    % reshape U to frame size
    U2 = reshape(U2, numel(ops.yrange), numel(ops.xrange), []);


    % write SVDs to disk
    if getOr(ops, {'writeSVDroi'}, 0)
        % project spatial masks onto raw data
        fid = fopen(ops.RegFile, 'r');
        ix = 0;
        Fs = zeros(ops.nSVDforROI, sum(ops.Nframes), 'single');
        while 1
            data = fread(fid,  Ly*Lx*nimgbatch, '*int16');
            if isempty(data)
                break;
            end
            data = single(data);
            data = reshape(data, Ly, Lx, []);

            % subtract off the mean of this batch
            data = bsxfun(@minus, data, mean(data,3));
            %     data = bsxfun(@minus, data, ops.mimg1);
            data = data(ops.yrange, ops.xrange, :);
            if ops.useGPU
                Fs(:, ix + (1:size(data,3))) = gpuBlockXtY(U, reshape(data, [], size(data,3)));
            else
                Fs(:, ix + (1:size(data,3))) = U' * reshape(data, [], size(data,3));
            end

            ix = ix + size(data,3);
        end
        fclose(fid);

        filePath = sm.getFileForPlane('svd_roi', ops.iplane);
        try
            save(filePath, 'U', 'Sv', 'Fs', 'ops', '-v6');
        catch
            save(filePath, 'U', 'Sv', 'Fs', 'ops');
        end
    end

end

