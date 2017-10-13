function dreg = register_F(data, ds)

[Lx, NN] = size(data);

Nx = ifftshift([-fix(Lx/2):ceil(Lx/2)-1]);

Nx = Nx' / Lx;

dph         = 2*pi*ds*Nx;

fdata       = fft(single(data));

dreg        = real(ifft(bsxfun(@times, fdata , exp(1i * dph))));



