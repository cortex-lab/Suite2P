dt = 1:1:120;

tsp = dcell{icell}.tspikes;
inds = repmat(tsp', numel(dt), 1) + repmat(dt', 1, numel(tsp));

F = dcell{icell}.F;
NT = numel(F);

nt = numel(dt);

basis = zeros(NT, nt);
bas0 = histc(tsp, .5:1:NT);
for k = 1:nt
   basis(k:end, k) = bas0(1:end-k+1);
end
%
basis(:,nt+1) = 1;

COV =basis'*basis;

coefs = COV\(basis' *F);
%
Frec = basis * coefs;

kernel = coefs(1:end-1);
%
plot(F)
hold all
plot(Frec)
hold off