function kernelS = realign_kernels(kernelS, tmax)

nt = size(kernelS, 1);

isub = linspace(1, nt, (nt-1)*10 + 1);

for j = 1:size(kernelS,2)
   ker = kernelS(:,j);
   kup = interp1(1:1:nt, ker, isub, 'spline');
   
   [~, imax] = max(kup);
   itmax = isub(imax); % the time of the max in units, needs to be at tmax
   
   
   kup = interp1(1:1:nt, ker, isub + itmax - tmax, 'spline', 0);
   
   kernelS(:,j) = kup(1:10:end);
end

