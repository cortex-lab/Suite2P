% compute bidirectional phase offsets from poor line scanning timing
function BiDiPhase = BiDiPhaseOffsets(data)

[Ly, Lx, nplanes, NT] = size(data);

% lines scanned one direction
yr1 = 2:2:floor(Ly/2)*2;
% lines scanned in other direction
yr2 = 1:2:floor(Ly/2)*2;

% compute phase correlation between lines in x-direction
eps0 = single(1e-6);
Nmax = min(50, NT);
d1 = fft(data(yr1,:,:,1:Nmax),[],2);
d2 = conj(fft(data(yr2,:,:,1:Nmax),[],2));
d1 = d1./(abs(d1) + eps0);
d2 = d2./(abs(d2) + eps0);

cc = ifft(d1 .* d2,[],2);
cc = fftshift(cc, 2);
cc = mean(mean(mean(cc,1),3),4);

% max shift of +/-5 pixels
[cx, ix] = max(cc(floor(Lx/2)+1 + [-5:5]));
ix       = ix - (6);

BiDiPhase = -1 * ix;



    