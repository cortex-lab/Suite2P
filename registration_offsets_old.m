function [ds, Corr] = registration_offsets_old(data, ops, remove_mean, varargin)

mimg = ops.mimg;

[Ly, Lx, NT] = size(data);

ys = repmat((1:Ly)', 1, Lx);
xs = repmat((1:Lx), Ly, 1);

ys = abs(ys - mean(ys(:)));
xs = abs(xs - mean(xs(:)));

mY      = max(ys(:)) - 4;
mX      = max(xs(:)) - 4;
slope   = 2;

Mask = 1./(1+exp((ys - mY)/slope)) ./(1+exp((xs - mX)/slope));

mMov        = mean(mimg(:));


xref = repmat(-2:+2, 5, 1);
yref = xref';

if ops.useGPU 
    data = gpuArray(data);
    Mask = gpuArray(single(Mask));
    mMov= gpuArray(single(mMov));
    mimg = gpuArray(single(mimg));
%     cc = gpuArray.zeros(size(data), 'single');
end
cc = zeros(size(data), 'single');
eps0 = single(1e-10);

cfft2mimg = conj(fft2(mimg));
cfft2mimgNorm = cfft2mimg./(eps0+abs(cfft2mimg));

for i = 1:NT
    J = mMov * (1-Mask) + Mask .* single(data(:,:,i));
    fdata =  fft2(J);
    if isfield(ops, 'PhaseCorrelation') && ~isempty(ops.PhaseCorrelation) ...
            && (ops.PhaseCorrelation==1)
        fdata =  fdata./(eps0+abs(fdata)); 
        f3 = fdata .* cfft2mimgNorm;
    else
        f3 = fdata .* cfft2mimg;
    end
    cc(:,:,i) =  gather(real(ifft2(f3)));
end

% at most lx, ly jitter
lx = 50; ly = 50;

xrange = [(Lx-lx+1):Lx 1:(lx+1)];
yrange = [(Ly-ly+1):Ly 1:(ly+1)];
cc =cc(yrange,xrange,:);

sig= 1;

ccs = my_conv(cc(:,:)', sig)';
ccs = reshape(ccs, size(cc));
ccs = permute(ccs, [2 1 3]);
ccs = my_conv(ccs(:,:)', sig)';
ccs = reshape(ccs, size(cc));
ccs = permute(ccs, [2 1 3]);

Corr = zeros(NT,1);
ds  = zeros(NT,2);
for i = 1:NT
     cc0 = cc(:,:,i);
     
    [dmax, iy] = max(ccs(:,:,i), [], 1);
    [mdmax, ix]    = max(dmax, [], 2);
    iy = iy(ix);
    
    if ops.SubPixel>1
        iy = min(max(iy, 3), 2*ly-1);
        ix = min(max(ix, 3), 2*lx-1);
        
        cczoom  = cc0(iy-2:iy+2, ix-2:ix+2);
        mcczoom = mean(mean(abs(cczoom)));
        ix = ix + mean(mean(xref .* cczoom))/mcczoom;
        iy = iy + mean(mean(yref .* cczoom))/mcczoom;
        Corr(i) = sum(sum(cczoom(:)));    
    else
        Corr(i) = mdmax;
    end
    ix = ix - lx - 1;
    iy = iy - ly - 1;
    
    if ops.SubPixel<Inf
        ix = round(ix * ops.SubPixel)/ops.SubPixel;
        iy = round(iy * ops.SubPixel)/ops.SubPixel;
    end
    
    ds(i,:) = [iy ix];
end

if sum(abs(ds(:))>200)>0
%     keyboard;
end

if remove_mean
   ds = ds - repmat(mean(ds,1), size(ds,1), 1); 
else
%     ds(Corr<ops.AlignNanThresh, :) = 0;
end




