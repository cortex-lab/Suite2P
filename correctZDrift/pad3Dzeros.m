% pad matrix in 3 dimensions
function B = pad3Dzeros(A,ny,nx,nz)

[Ny, Nx, Nz] = size(A);
pD = [floor(ny/2) floor(nx/2) floor(nz/2)];
pOff = [ny nx nz] - pD*2;

B = padarray(A,pD+pOff,0,'pre');
B = padarray(B,pD,0,'post');


end