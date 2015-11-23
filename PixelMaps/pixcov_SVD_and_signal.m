function pixCov = pixcov_SVD_and_signal(U, V, R)

R = R(:);

R = R - mean(R);
V = V - repmat(mean(V,2), 1, size(V,2));

NT   = size(V,2);
svdC = (V * R)/NT;

[Ly, Lx, nSVD] = size(U);

pixCov = reshape(reshape(U, Ly*Lx, nSVD)  * svdC, Ly, Lx);

end