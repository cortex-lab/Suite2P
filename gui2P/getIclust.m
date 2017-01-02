function [iclust, lambda, lambda0] = getIclust(stat, cl)

iclust = zeros(cl.Ly, cl.Lx);
lambda = zeros(cl.Ly, cl.Lx);
lambda0 = zeros(cl.Ly, cl.Lx);

for j = 1:numel(stat)
   ipix = stat(j).ipix;
   
   inew = stat(j).lambda(:)>lambda(ipix) + 1e-6;
%    lambda(ipix(inew)) = stat(j).lam(inew);
%    lambda0(ipix(inew)) = stat(j).lambda(inew);
   
   lambda(ipix(inew)) = stat(j).lambda(inew);
   lambda0(ipix(inew)) = stat(j).lam(inew);
   
   iclust(ipix(inew)) = j;    
end
