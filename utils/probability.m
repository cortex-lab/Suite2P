
function Ypred = probability(train_stats, train_iscell, test_stats)

nodes = 100;
[nroi, nstats] = size(train_stats);
[ssort, isort] = sort(train_stats, 1);
ix = round(linspace(1, nroi,nodes));
grid = ssort(ix,:);
grid(end,:) = Inf;
p = zeros(nodes-1, nstats);
for j = 1:nodes-1
   for k = 1:nstats
      p(j,k) = mean(train_iscell(isort(ix(j):(ix(j+1)-1), k))); 
   end
end
p = my_conv2(p, 2, 1);
logp = get_logp(test_stats, grid, p);
logp = sum(logp,2);
Ypred = 1./(1+exp(-logp));

end

function logp = get_logp(test_stats, grid, p)
[nroi, nstats] = size(test_stats);
logp = zeros(nroi, nstats);
for n = 1:nstats
   x = test_stats(:,n);
   x(x<grid(1,n)) = grid(1,n);
   x(x>grid(end,n)) = grid(end,n);
   [~, ibin] = histc(x, grid(:,n));
   logp(:,n) = log(p(ibin, n) + 1e-6) - log(1-p(ibin, n) + 1e-6);
end

end