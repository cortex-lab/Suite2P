function [zpos, FcellZ, FcellNeuZ] = writeZPostoBin(ops, PtoZ, MimgZ)

% use registered binary files from suite2P
% make sure that ops.DeleteBin = 0

zpos   = zeros(sum(ops.Nframes), ops.nplanes,'single');
%%
ipl = [];
for j = 1:length(PtoZ)
    if ~isempty(PtoZ{j})
        ipl = [ipl j];
    end
end

for iplane = ipl
    
    planefile = sprintf('%s/F_%s_%s_plane%d.mat', ops.ResultsSavePath, ...
        ops.mouse_name, ops.date, iplane);
    dat = load(planefile);
    
    %%
    Taff = PtoZ{iplane};
    
    % cut through z-stack where mean plane should be +/- zspread
    Ny    =  numel(dat.ops.yrange);
    Nx    =  numel(dat.ops.xrange);
    
    mimg  = dat.ops.mimg1(dat.ops.yrange, dat.ops.xrange);
    
    zspread  = 10;
    Zaligned = cutZstack(MimgZ, Taff, Ny, Nx, zspread);
    
    % rescale z-stack
    Zaligned   = single(Zaligned) / mean(single(Zaligned(:))) * mean(single(mimg(:)));
    
    % check that cut is centered on mean plane
    m1 = fft(fft(Zaligned,[],1),[],2);
    eps0 = single(1e-20);
    m1 = m1./(abs(m1)+eps0);
    m2 = fft(fft(mimg,[],1),[],2);
    m2 = m2./(abs(m2)+eps0);
    if ops.useGPU
        m1 = gpuArray(single(m1));
        m2 = gpuArray(single(m2));
    end
    
    [cx, ix, cZ] = ZRegPlane(m1,m2,[1:size(m1,3)],ops.useGPU);
    clear m1 m2;
    ix(3)      = zspread + 1 - ix(3);
    fprintf('PLANE %d \n x offset: %2.2f; y offset: %2.2f; z offset: %2.2f\n',iplane, ix(1), ix(2), ix(3));
    
    
    subpixel     = 50;
    iZ           = [-zspread : 1/subpixel : zspread];
    
    % upsample z-stack
    sigL         = 1; % kernel width in pixels
    zt           = [-zspread : zspread];
    Kx           = kernelD(zt,zt,sigL);
    Kg           = kernelD(iZ,zt,sigL);
    Kmat         = single(Kg/Kx);
    Kmat         = bsxfun(@rdivide, Kmat, sum(Kmat,2));
    Zp           = reshape(permute(Zaligned, [3 1 2]), numel(zt), []);
    Zups         = (Kmat * Zp)';
    %Zups         = bsxfun(@plus, (Kmat * bsxfun(@minus, Zp, mean(Zp,1))), mean(Zp,1))';
    Zups         = reshape(Zups, Ny, Nx, numel(iZ));
    %Zwhite = reshape(Zups, [], size(Zups,3));
    
    % whiten z-stack (take ifft of phase of fft)
    m1 = fft(fft(Zaligned,[],1),[],2);
    eps0 = single(1e-20);
    m1 = m1./(abs(m1)+eps0);
    Zwhite = real(ifft(ifft(m1, [], 1), [], 2));
    Zwhite = reshape(Zwhite, [], size(Zwhite,3));
    Zwhite = gpuArray(single(Zwhite));
    
    
    % check mean image
    clf;
    subplot(1,2,1),
    imagesc(mimg,[0 8000]);
    title(sprintf('plane %d',iplane));
    subplot(1,2,2),
    imagesc(Zaligned(:,:,zspread+1),[0 8000]);
    title('stretched z-stack');
    drawnow;
    
    %%
    for j = find(nPos>0)'
        clf;
        colormap('gray');
        subplot(1,2,1),
        imagesc(Zups(:,:,j),[500 5000]);
        title(j)
        subplot(1,2,2),
        imagesc(mimgPos(:,:,j),[500 5000]);
        title(nPos(j));
        drawnow
        pause;
    end
    
    
    %% open registered binary file and correlate with z-stack to find z shifts
    ops0 = dat.ops;
    ops0.Ly = Ny;
    ops0.Lx = Nx;
    ops0.numBlocks = [3 3];
    ops0.blockFrac = 0.4;
    ops0 = MakeBlocks(ops0);
    yBL = ops0.yBL;
    xBL = ops0.xBL;
    clear iBL;
    for j = 1:length(yBL)
        iBL{j} = sub2ind([Ny Nx],  repmat(yBL{j}',1,numel(xBL{j})), repmat(xBL{j},numel(yBL{j}),1));
    end
    
    xyMask = ops0.xyMask;
    
    %%
    fid = fopen(dat.ops.RegFile, 'r');
    
    Ly  = dat.ops.Ly;
    Lx  = dat.ops.Lx;
    NT  = size(dat.ops.DS,1);
    Nbatch = 250 / round(Ly/512);
    ix0    = 0;
    tic;
    kk = 0;
        
    % array for z-position of planes
    mimgPos = single(Zups);
    nPos    = ones(numel(iZ), 1, 'single');
    while ix0 < NT
        indxr = ix0 + (1:Nbatch);
        ix0   = ix0 + Nbatch;
        indxr(indxr > NT) = [];
        
        data  = fread(fid,  Ly*Lx*length(indxr), '*int16');
        data  = reshape(data, Ly, Lx, []);
        data  = data(dat.ops.yrange, dat.ops.xrange,:);
        data0 = single(data);
        
        if ops.useGPU
            data = gpuArray(single(data));
        end
        
        % whiten data
        m2     = fft(fft(data, [], 1), [], 2);
        data   = real(ifft(ifft(m2, [], 1), [], 2));
        data   = reshape(data, [], size(data,3));
        
        % correlate with z-stack
        %cc     = Zwhite' * data;
        %if ops.useGPU
        %    cc = gather(cc);
        %end
                
        % try block correlations
        cc     = zeros(size(Zwhite,2), size(data,2), length(iBL));
        for j = 1:length(iBL)
            cc0 = Zwhite(iBL{j},:)' * data(iBL{j},:);
            if ops.useGPU
                cc0 = gather(cc0);
            end
            cc(:,:,j) = cc0;
        end
        
        % interpolate and find max position in z-stack
        ccz0          = interp1([-zspread:zspread]', cc, iZ(:), 'spline');
        [cmax,izmax0] = max(ccz0);
        izmax0        = squeeze(round(izmax0));
        [~,imax]      = max(squeeze(cmax),[],2);
        
        for j = 1:length(indxr)
            % which block has max correlation? constrain others to be +/-3
            ccz       = squeeze(ccz0(:,j,:));
            ij        = izmax0(j,imax(j));
            cinds     = max(1, ij - subpixel*2) : min(numel(iZ), ij + subpixel*2);  
            ccz       = ccz(cinds,:);
            [~,izmax] = max(ccz,[],1);
            izmax     = izmax + cinds(1) - 1;
            
            zpos(indxr(j), :, iplane)  = iZ(izmax);
            zind(indxr(j), :, iplane)  = izmax;
            
            izmax = round(mean(izmax,2));
            mimgPos(:, :, izmax) = mimgPos(:, :, izmax) + data0(:, :, j);
            nPos(izmax)          = nPos(izmax) + 1;
        end
        
        if rem(ix0, 6000)==0
            fprintf('Frame %d done in time %2.2f \n', ix0, toc)
        end
    end
    fclose(fid);
    
    mimgPos = bsxfun(@rdivide, mimgPos, permute(nPos, [2 3 1]));
    
    %% take z-shifts and make FcellZ and FcellNeu
    Nk       = numel(dat.stat); % all ROIs
    ops      = dat.ops;
    
    stat = getNonOverlapROIs(stat, Ny, Nx);
    % create cell masks and cell exclusion areas
    [stat, cellPix, cellMasks] = createCellMasks(dat.stat, Ny, Nx);
    % create surround neuropil masks
    [ops, neuropMasks] = createNeuropilMasks(ops, stat, cellPix);
    % add surround neuropil masks to stat
    for k = 1:Nk
        stat(k).ipix_neuropil = find(squeeze(neuropMasks(k,:,:))>0);
    end
    % convert masks to sparse matrices for fast multiplication
    neuropMasks = sparse(double(neuropMasks(:,:)));
    cellMasks   = sparse(double(cellMasks(:,:)));
    
    %%
    ix0 = 0;
    Nbatch = 2000;
    F      = zeros(numel(stat), NT, 'single');
    Fneu      = zeros(numel(stat), NT, 'single');
    mimgPos   = reshape(double(mimgPos), [], size(mimgPos,3));
    while ix0 < NT
        indxr = ix0 + (1:Nbatch);
        ix0   = ix0 + Nbatch;
        indxr(indxr > NT) = [];
        zindxr            = round(squeeze(mean(zind(indxr,:,iplane),2)));
        Ftemp             = cellMasks * mimgPos(:,zindxr);
        F(:, indxr)       = single(Ftemp);
        Ftemp             = neuropMasks * mimgPos(:,zindxr);
        Fneu(:,indxr)     = single(Ftemp);
    end
    mimgPos    = reshape(single(mimgPos), Ny, Nx, []);
    %%
    
    
    csumNframes = [0 cumsum(ops.Nframes)];
    FcellZ       = cell(1, length(ops.Nframes));
    FcellNeuZ    = cell(1, length(ops.Nframes));
    for i = 1:length(ops.Nframes)
        FcellZ{i}     = F(:, csumNframes(i) + (1:ops.Nframes(i)));
        FcellNeuZ{i}  = Fneu(:, csumNframes(i) + (1:ops.Nframes(i)));
    end

    
    %%
    fz  = fopen(sprintf('%s_Z.bin',dat.ops.RegFile(1:end-4)), 'w');
    ix0    = 0;
    tic;
    while ix0 < NT
        indxr = ix0 + (1:Nbatch);
        ix0   = ix0 + Nbatch;
        indxr(indxr > NT) = [];
        
        zd  = reshape(mimgPos(:,:,zind(indxr,:,iplane)),Ny*Nx,numel(indxr),[]);
        zd  = bsxfun(@times, permute(zd, [1 3 2]), xyMask);
        zd  = reshape(squeeze(sum(zd, 2)), Ny, Nx, numel(indxr));
        
        zdata        = zeros(Ly, Lx, length(indxr), 'int16');
        zdata(dat.ops.yrange, dat.ops.xrange, :) = ...
            int16(round(zd));
        fwrite(fz, zdata, 'int16');            
        
        if rem(ix0, 6000)==0
            fprintf('Frame %d done in time %2.2f \n', ix0, toc)
        end
    end
    fclose(fz);
    
end
