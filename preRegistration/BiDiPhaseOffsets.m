function BiDiPhase = BiDiPhaseOffsets(data)

[Ly, Lx, nplanes, NT] = size(data);

yrange = 2:2:Ly;

%
yr1 = 2:2:Ly;
yr2 = 1:2:Ly;

eps0 = single(1e-20);
d1 = fft(data(yr1,:,:,:),[],2);
d2 = fft(data(yr2,:,:,:),[],2);
d1 = d1./(abs(d1) + eps0);
d2 = d2./(abs(d2) + eps0);

cc = ifft(d1 .* conj(d2),[],2);
cc = fftshift(cc, 2);
cc = mean(mean(mean(cc,1),3),4);

[cx, ix] = max(cc);
ix       = ix - (floor(Lx/2) + 1);

BiDiPhase = ix;



    