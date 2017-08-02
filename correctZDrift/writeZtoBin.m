function zpos = writeZtoBin(ops, PtoZ, MimgZ)

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
    ly    =  numel(dat.ops.yrange);
    lx    =  numel(dat.ops.xrange);
    
    mimg  = dat.ops.mimg1(dat.ops.yrange, dat.ops.xrange);
    
    zspread  = 10;
    Zaligned = cutZstack(MimgZ, Taff, ly, lx, zspread);
    
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
    
    %% kriging interpolation so that Zaligned is subpixel for shifts
    % compute kernels for regression
    subpixel     = 50;
    iZ           = [-zspread : 1/subpixel : zspread];
    sigL         = 1; % kernel width in pixels
    zt           = [-zspread : zspread];
    Kx           = kernelD(zt,zt,sigL);
    Kg           = kernelD(iZ,zt,sigL);
    Kmat         = single(Kg/Kx);
    Kmat         = bsxfun(@rdivide, Kmat, sum(Kmat,2));
    Zp           = reshape(permute(Zaligned, [3 1 2]), numel(zt), []);
    Zups         = (Kmat * Zp)';
    %Zups         = bsxfun(@plus, (Kmat * bsxfun(@minus, Zp, mean(Zp,1))), mean(Zp,1))';
    Zups         = reshape(Zups, ly, lx, numel(iZ));
    %Zwhite = reshape(Zups, [], size(Zups,3));
    
    %%
    
    m1 = fft(fft(Zaligned,[],1),[],2);
    eps0 = single(1e-20);
    m1 = m1./(abs(m1)+eps0);
    % whiten z-stack (take ifft of phase of fft)
    Zwhite = real(ifft(ifft(m1, [], 1), [], 2));
    Zwhite = reshape(Zwhite, [], size(Zwhite,3));
    Zwhite = gpuArray(single(Zwhite));
    
    %%
    %
    %     %icell = ceil(rand*1100);
    %     Z2 = reshape(Zups,[],size(Zups,3));
    %     imask = stat(icell).lam/sum(stat(icell).lam);
    %     cz = Z2(stat(icell).ipix,:)'*imask;
    %
    %     figure(1)
    %     subplot(1,3,1),
    %     hist(cz(floor(zpos(:,3)*subpixel)+subpixel*zspread+1),40);
    %     subplot(1,3,2),
    %     plot(cz);
    %     subplot(1,3,3),
    %     plot(floor(zpos(:,3)*subpixel)+subpixel*zspread+1);
    %
    %     figure(2)
    %     plot(cz(floor(zpos(:,3)*subpixel)+subpixel*zspread+1))
    %     %%
    %     for j = 1:size(Zups,3)
    %         clf;
    %         subplot(1,2,1),
    %         imagesc(Zups(:,:,j),[0 5000]);
    %         title(j);
    %         subplot(1,2,2),
    %         imagesc(Zaligned(:,:,floor(j/subpixel)+1),[0 5000]);
    %         drawnow;
    %         pause(.025);
    %     end
    %     %%
    
    % check mean image
    clf;
    subplot(1,2,1),
    imagesc(mimg,[0 8000]);
    title(sprintf('plane %d',iplane));
    subplot(1,2,2),
    imagesc(Zaligned(:,:,zspread+1),[0 8000]);
    title('stretched z-stack');
    drawnow;
    
    %% open registered binary file and correlate with z-stack to find z shifts
    
    
    fid = fopen(dat.ops.RegFile, 'r');
    
    fz  = fopen(sprintf('%s_Z.bin',dat.ops.RegFile(1:end-4)), 'w');
    
    Ly  = dat.ops.Ly;
    Lx  = dat.ops.Lx;
    NT  = size(dat.ops.DS,1);
    Nbatch = 250 / round(Ly/512);
    ix0    = 0;
    tic;
    kk = 0;
    zex = zeros(Ly, Lx, 1);
    while ix0 < NT
        indxr = ix0 + (1:Nbatch);
        ix0   = ix0 + Nbatch;
        indxr(indxr > NT) = [];
        
        data  = fread(fid,  Ly*Lx*length(indxr), '*int16');
        data  = reshape(data, Ly, Lx, []);
        data  = data(dat.ops.yrange, dat.ops.xrange,:);
        
        if ops.useGPU
            data = gpuArray(single(data));
        end
        
        % whiten data
        m2     = fft(fft(data, [], 1), [], 2);
        data   = real(ifft(ifft(m2, [], 1), [], 2));
        data   = reshape(data, [], size(data,3));
        
        % correlate with z-stack
        cc     = Zwhite' * data;
        if ops.useGPU
            cc = gather(cc);
        end
        
        % interpolate and find max position in z-stack
        ccz          = interp1([-zspread:zspread]', cc, iZ(:), 'spline');
        [cmax,izmax] = max(ccz);
        zpos(indxr,iplane)  = iZ(izmax);
        zind         = izmax;
        
        
        % put zdata into 512 x 512 format
        zdata        = zeros(Ly, Lx, length(indxr), 'int16');
        for j = 1:length(indxr)
            zdata(dat.ops.yrange, dat.ops.xrange, j) = ...
                int16(round(Zups(:,:,zind(j))));
        end
        fwrite(fz,  zdata, 'int16');
        
        if rem(ix0, 1000) == 0
            kk=kk+1;
            zex(:,:,kk) = zdata(:,:,j);
        end
        
        if rem(ix0, 6000)==0
            fprintf('Frame %d done in time %2.2f \n', ix0, toc)
        end
    end
    fclose(fid);
    fclose(fz);
end
