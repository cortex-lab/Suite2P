function [bestphi, amps, widths, A] = get_vonMisesFit(R)
% assumes orientations showed are linearly space from 0 to 180, excluding
% 180. Discretizes results, for example orientations in 3 degree steps. 

[~, Nori] = size(R); 

oritest = repmat(linspace(0, pi, 12*5+1), 24, 1);
oritest = oritest(:,1:12*5);
ktest = repmat(linspace(0, 4, 24)', 1, 12*5);

oribase = linspace(0, pi, Nori+1);
oribase = oribase(1:end-1);

f = zeros([size(oritest), length(oribase)]);
for i = 1:length(oribase)
    f(:,:,i) = exp(ktest.*cos(2*(oribase(i) - oritest)))./(2*pi*besseli(0, ktest));
end

f = reshape(f, [], length(oribase));
fnorm = sum(f.^2,2).^.5;
f = f./repmat(fnorm, 1, length(oribase));

coefs = f * R'; 

[M, imax] = max(coefs,[], 1);

bestphi = oritest(imax);
widths  = 1./ktest(imax).^.5;
amps = M'.*fnorm(imax); 


A = abs(rem(repmat(bestphi, Nori, 1) - repmat(oribase', 1, size(R,1)), pi));
A = min(A, pi-A);
[~, A] = min(A, [], 1);