function [zpos] = zposPlane(ops, Taff, MimgZ, subpixel)

% use Z-stack and registered binary files from suite2P to compute Z-position of plane
% make sure that ops.DeleteBin = 0

zpos   = zeros(sum(ops.Nframes), 1,'single');
    
% cut through z-stack where mean plane should be +/- zspread
Ny    =  numel(ops.yrange);
Nx    =  numel(ops.xrange);

mimg  = ops.mimg1(ops.yrange, ops.xrange);

%%
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
fprintf('>>> x offset: %2.2f; y offset: %2.2f; z offset: %2.2f\n', ix(1), ix(2), ix(3));

% whiten z-stack (take ifft of phase of fft)
m1 = fft(fft(Zaligned,[],1),[],2);
eps0 = single(1e-20);
m1 = m1./(abs(m1)+eps0);
Zwhite = real(ifft(ifft(m1, [], 1), [], 2));
clear m1;
Zwhite = reshape(Zwhite, [], size(Zwhite,3));
if ops.useGPU
    Zwhite = gpuArray(single(Zwhite));
end
        
% check mean image
clf;
subplot(1,2,1),
imagesc(mimg,[0 8000]);
subplot(1,2,2),
imagesc(Zaligned(:,:,zspread+1),[0 8000]);
title('stretched z-stack');
drawnow;
    
%%
iZ       = [-zspread : 1/subpixel : zspread];
fid = fopen(ops.RegFile, 'r');
    
Ly  = ops.Ly;
Lx  = ops.Lx;
NT  = sum(ops.Nframes);
Nbatch = 250 / round(Ly/512);
ix0    = 0;
tic;
        
% array for z-position of planes
nPos    = ones(numel(iZ), 1, 'single');
while ix0 < NT
    indxr = ix0 + (1:Nbatch);
    ix0   = ix0 + Nbatch;
    indxr(indxr > NT) = [];
    
    data  = fread(fid,  Ly*Lx*length(indxr), '*int16');
    data  = reshape(data, Ly, Lx, []);
    data  = data(ops.yrange, ops.xrange, :);
    if ops.useGPU
        data = gpuArray(single(data));
    else
        data = single(data);
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
    ccz           = interp1([-zspread:zspread]', cc, iZ(:), 'spline');
    [cmax,izmax0] = max(ccz);
    izmax         = squeeze(round(izmax0));
    zpos(indxr)  = iZ(izmax);
        
    if rem(ix0, 6000)==0
        fprintf('Frame %d done in time %2.2f \n', ix0, toc)
    end
end
fclose(fid);
