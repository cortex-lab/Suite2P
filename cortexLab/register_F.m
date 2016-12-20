function dreg = register_F(data, ds)

[Lx, NN] = size(data);

if numel(ds)==1
   ds = repmat(ds, NN,1); 
end
Nx = ifftshift([-fix(Lx/2):ceil(Lx/2)-1]);

Nx = Nx' / Lx;

dreg  = zeros(size(data), 'like', data);

for i = 1:NN
    dph         = 2*pi*ds(i)*Nx;
    fdata       = fft(single(data(:,i)));
    dreg(:,i) = real(ifft(fdata .* exp(1i * dph)));
end

