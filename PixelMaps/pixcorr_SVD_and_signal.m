function pixCorr = pixcorr_SVD_and_signal(U, V, R)

% Computes pixel-level correlation with secondary signal R, such as a
% running trace or a pupil trace. U and V are the spatial and temporal components
% of the SVD respectively, such that the raw fluorescence is F = U * V. 
% This script produces an array of correlation coefficients
% pixCorr of size Ly by Lx. The number of SVD components used is whatever
% size(U,3) is. 

[Ly, Lx, nSVD] = size(U);

if nSVD>size(V,1)
    error('the number of dimensions in V is smaller than the number of dimensions in U');
end
V = V(1:nSVD, :);

pixCov = pixcov_SVD_and_signal(U, V, R); % covariance between each pixel and trace R
pixVar = pixvar_SVD(U, V);  % variance of each pixel, computed from SVDs
varR   = var(R, 1); % variance of R, normalized by full number of samples

pixCorr = pixCov ./(varR * pixVar).^.5;

end