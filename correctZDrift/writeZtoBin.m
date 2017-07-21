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
    
    Taff = PtoZ{iplane};
    
    % cut through z-stack where mean plane should be +/- zspread
    ly    =  numel(dat.ops.yrange);
    lx    =  numel(dat.ops.xrange);
    
    mimg  = dat.ops.mimg1(dat.ops.yrange, dat.ops.xrange);
    
    zspread  = 10;
    Zaligned = cutZstack(MimgZ, Taff, ly, lx, zspread);
    
    % rescale z-stack
    Zaligned   = Zaligned / mean(Zaligned(:)) * mean(mimg(:));
        
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
    ix(3)      = zspread + 1 - ix(3); 
    fprintf('PLANE %d \n x offset: %2.2f; y offset: %2.2f; z offset: %2.2f\n',iplane, ix(1), ix(2), ix(3));
    
    
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
    
    % whiten z-stack (take ifft of phase of fft)
    Zwhite = real(ifft(ifft(m1, [], 1), [], 2));
    Zwhite = reshape(Zwhite, [], size(Zwhite,3));
   
    fid = fopen(dat.ops.RegFile, 'r');
    
    fz  = fopen(sprintf('%s_Z.bin',dat.ops.RegFile(1:end-4)), 'w');
    
    Ly  = dat.ops.Ly;
    Lx  = dat.ops.Lx;
    NT  = size(dat.ops.DS,1);
    Nbatch = 250 / round(Ly/512);
    ix0    = 0;
    tic;
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
        dwhite = real(ifft(ifft(m2, [], 1), [], 2));
        dwhite = reshape(dwhite, [], size(dwhite,3));
        
        % correlate with z-stack
        cc     = Zwhite' * dwhite;
        if ops.useGPU
            cc = gather(cc);
        end
        % interpolate and find max position in z-stack
        iZ           = [-zspread: .01 : zspread];
        cZi          = interp1([-zspread:zspread], cc, iZ, 'spline');
        [cmax,izmax] = max(cZi);
        zpos(indxr,iplane)  = iZ(izmax);
        zind         = round(iZ(izmax)) + zspread + 1;
        
        % put zdata into 512 x 512 format
        zdata        = zeros(Ly, Lx, length(indxr), 'int16');
        for j = 1:length(indxr)
            zdata(dat.ops.yrange, dat.ops.xrange, j) = ...
                int16(round(Zaligned(:,:,zind(j))));
        end
        fwrite(fz,  zdata, 'int16');
        
        if rem(ix0, 6000)==0
            fprintf('Frame %d done in time %2.2f \n', ix0, toc)
        end
    end
    fclose(fid);
    fclose(fz);
end
