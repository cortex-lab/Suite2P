function [iclust, lambda] = getIclust(stat, cl)

iclust = zeros(cl.Ly, cl.Lx);
lambda = zeros(cl.Ly, cl.Lx);

for j = 1:numel(stat)
   ipix = stat(j).ipix;
   
   inew = stat(j).lambda(:)>lambda(ipix) + 1e-6;
   lambda(ipix(inew)) = stat(j).lambda(inew);
   
   iclust(ipix(inew)) = j;    
end
