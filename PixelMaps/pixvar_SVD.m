function pixVar = pixvar_SVD(U, V)

V = V - repmat(mean(V,2), 1, size(V,2));

NT   = size(V,2);
svdV = (V * V')/NT;

[Ly, Lx, nSVD] = size(U);

pixVar = reshape(U, Ly*Lx,nSVD) * svdV; 

pixVar = sum(pixVar .* reshape(U, Ly*Lx,nSVD), 2);
pixVar = reshape(pixVar, Ly, Lx);

end