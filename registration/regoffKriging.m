function [dv,corr,regdata] = regoffKriging(data, ops, removeMean)

refImg = ops.mimg;
subpixel = getOr(ops, {'subPixel' 'SubPixel'}, 10); % subpixel factor
useGPU = getOr(ops, 'useGPU', false);
phaseCorrelation = getOr(ops, {'phaseCorrelation' 'PhaseCorrelation'}, true);
maxregshift = getOr(ops, 'maxregshift', 50);

[ly lx nFrames] = size(data);

% return registered movie
if nargout > 2 % translation required
    translate = true;
    fy = ifftshift((-fix(ly/2):ceil(ly/2) - 1)/ly)';% freq along first dimension
    fx = ifftshift((-fix(lx/2):ceil(lx/2) - 1)/lx); % freq along second dimension
else
    translate = false;
end

eps0          = single(1e-20);
if useGPU
    batchSize = getBatchSize(ly*lx);
    eps0      = gpuArray(eps0);
    refImg    = gpuArray(single(refImg));
else
    batchSize = 1000;
end

% fft of reference image
m2 = fftn(refImg);
if phaseCorrelation
    m2 = m2./(eps0 + abs(m2));
end

% allow max shifts +/- lcorr
lpad   = 3;
lcorr  = min(maxregshift, floor(min(ly,lx)/2)-lpad);

% only need a small kernel +/- lpad for smoothing
[x1,x2] = ndgrid([-lpad:lpad]);
xt = [x1(:) x2(:)]';
if useGPU
    xt = gpuArray(single(xt));
end
    
if ops.kriging
    % compute kernels for regression
    sigL     = 0.75; % kernel width in pixels
    Kx = kernelD(xt,xt,sigL*[1;1]);
    linds = [-lpad:1/subpixel:lpad];
    [x1,x2] = ndgrid(linds);
    xg = [x1(:) x2(:)]';
    if useGPU
        xg = gpuArray(single(xg));
    end
    Kg = kernelD(xg,xt,sigL*[1;1]);
    Kmat = Kg/Kx;
end

% loop over batches
dv = zeros(nFrames, 2);
corr = zeros(nFrames, 1);
if translate
    regdata = zeros(ly, lx, nFrames, 'single');
end

nBatches = ceil(nFrames/batchSize);
for bi = 1:nBatches
    fi = (bi - 1)*batchSize + 1:min(bi*batchSize, nFrames);
    
    if useGPU
        batchData = gpuArray(single(data(:,:,fi)));
    else
        batchData = single(data(:,:,fi));
    end
    
    m1 = fft(fft(batchData,[],1),[],2);
    if phaseCorrelation
        m1 = m1./(eps0 + abs(m1));
    end
    
    % compute correlation matrix
    cc = real(ifft(ifft(bsxfun(@times,m1,conj(m2)),[],1),[],2));
    cc = fftshift(fftshift(cc,1),2);
        
    %% subpixel registration
    if subpixel > 1
        % kriging subpixel
        for j = 1:numel(fi)
            % allow only +/- lcorr shifts
            cc0     = cc(floor(ly/2)+1+[-lcorr:lcorr],floor(lx/2)+1+[-lcorr:lcorr],j);
            [~,ix] = max(cc0(:));
            [ix1,ix2] = ind2sub(size(cc0),ix);
            % max point within +/- lcorr wrt whole matrix
            mxpt    = [ix1+floor(ly/2) ix2+floor(lx/2)]-lcorr;
            
            % matrix +/- lpad surrounding max point
            ypad    = max(1, mxpt(1)-lpad) : min(ly, mxpt(1)+lpad);
            iy0     = min(0, mxpt(1)-lpad-1);
            xpad    = max(1, mxpt(2)-lpad) : min(lx, mxpt(2)+lpad);
            ix0     = min(0, mxpt(2)-lpad-1);
            if useGPU
                ccmat   = gpuArray.zeros(2*lpad+1,2*lpad+1,'single');
            else
                ccmat   = zeros(2*lpad+1,2*lpad+1,'single');
            end
            ccmat(iy0+[1:numel(ypad)], ix0+[1:numel(xpad)])   = cc(ypad,xpad,j);
            
            if ops.kriging
                % regress onto subsampled grid
                ccb    = Kmat * ccmat(:);
                
                % find max of grid
                [cx,ix] = max(ccb(:));
                [ix11,ix21] = ind2sub(numel(linds)*[1 1],ix);
                mdpt    = floor(numel(linds)/2)+1;
                dv0   = ([ix11 ix21] - mdpt)/subpixel + mxpt - [floor(ly/2) floor(lx/2)] - 1;
            else
                yshift = xt(1,:) * ccmat(:);
                xshift = xt(2,:) * ccmat(:);
                dv0    = [yshift xshift]/sum(ccmat(:)) + mxpt - [floor(ly/2) floor(lx/2)] - 1;
                if isfinite(subpixel)
                    dv0 = round(dv0 * subpixel) / subpixel;
                end
                cx     = max(ccmat(:));
            end
            dv(fi(j),:) = gather_try(dv0);
            corr(fi(j))  = gather_try(cx);
            
        end
        
    % otherwise just take peak of matrix
    else
        cc0     = cc(floor(ly/2)+1+[-lcorr:lcorr],floor(lx/2)+1+[-lcorr:lcorr],:);
        [cmax,iy]  = max(cc0,[],1);
        [cx, ix]   = max(cmax,[],2);
        iy = reshape(iy(sub2ind([size(iy,2) size(iy,3)], ix(:), (1:size(iy,3))')),...
            1, 1, []);
        
        dv0 = [iy(:)-lcorr ix(:)-lcorr]-1;
        dv(fi,:)  = gather_try(dv0);
        corr(fi) = gather_try(cx(:));
    end
    
    % shift movie
    if translate
        phaseShift = bsxfun(@times,...
            exp(1j*2*pi*bsxfun(@times, fy, dv0(:,1))),... y rotation
            exp(1j*2*pi*bsxfun(@times, fx, dv0(:,2)))); % x rotation
        res = real(ifft2(fft2(batchData).*phaseShift));
        regdata(:,:,fi) = gather_try(res);
    end
    
end



