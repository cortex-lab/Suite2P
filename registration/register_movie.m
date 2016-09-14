function dreg = register_movie(data, ops, ds)

orig_class = class(data);

if ops.useGPU
    data = gpuArray(single(data));
end
[Ly, Lx, NT] = size(data);

Ny = ifftshift([-fix(Ly/2):ceil(Ly/2)-1]);
Nx = ifftshift([-fix(Lx/2):ceil(Lx/2)-1]);
[Nx,Ny] = meshgrid(Nx,Ny);
Nx = Nx / Lx;
Ny = Ny / Ly;

if ops.useGPU
    dreg = gpuArray.zeros(size(data), orig_class);
else
    dreg = zeros(size(data), orig_class);
end


if ops.useGPU
    ds = gpuArray(ds);
    Nx = gpuArray(single(Nx));
    Ny = gpuArray(single(Ny));
end

for i = 1:NT
    dph         = 2*pi*(ds(i,1)*Ny + ds(i,2)*Nx);
    fdata       = fft2(single(data(:,:,i)));
    dreg(:,:,i) = real(ifft2(fdata .* exp(1i * dph)));
end

if ops.useGPU
    dreg = gather_try(dreg);
end
