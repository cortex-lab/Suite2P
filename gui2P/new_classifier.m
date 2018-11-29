function h = new_classifier(h)

[~, filename1] = fileparts(h.dat.cl.fpath);
set(h.text50,'String', [filename1 '.mat']);

% this load iscell, st and statLabels
load(h.dat.cl.fpath);
stat = h.dat.stat;

st0 = zeros( numel(stat), numel(statLabels));
for j = 1:numel(statLabels)
    if isfield(stat, statLabels{j})        
        st0(:,j) = eval(sprintf('[stat.%s]', statLabels{j}));
    else
        error('could not find field %s in stat', statLabels{j})
    end
end
st0(isnan(st0)) = 2;

Ypred = probability(st, icell, st0);

for j = 1:length(stat)
	if h.init==0
		stat(j).iscell = Ypred(j) > h.dat.cl.threshold;
	end
   stat(j).cellProb  = Ypred(j); 
end

h.dat.stat      = stat;
h.st0           = st0;
h.icell           = icell;
h.statLabels    = statLabels;

end

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