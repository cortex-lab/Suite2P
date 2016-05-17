%%%% shifts a single point in 2D fft and gives new coordinates
%%%% shiftType = 1 is fft, else inverse
function ixS = FftShift2D(ix,nX,nY,shiftType)

ix0 = ix;
ix = round(ix);

ixcc       = zeros(nX,nY);
ixcc(ix(1),ix(2)) = 1;
if shiftType == 1
    ixcc       = fftshift(ixcc);
    ixcc0      = fftshift(ix-round(ix0));
else
    ixcc       = ifftshift(ixcc);
    ixcc0      = ifftshift(ix-round(ix0));
end

ixS         = [0 0];
[ixS(1), ixS(2)] = find(ixcc==1);

ixS = ixS + ixcc0;