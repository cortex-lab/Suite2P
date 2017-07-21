function [stat, F, Fneu] = getActivityStats(ops, stat, F, Fneu)

indNoNaN    = find(~ops.badframes);
ix          = cumsum(~ops.badframes) + 1;
ix          = ix(ops.badframes);
ix(ix>numel(indNoNaN))  = numel(indNoNaN);

F(:, ops.badframes)  = F(:,    indNoNaN(ix));
Fneu(:, ops.badframes)  = Fneu(:, indNoNaN(ix));


coefNeu = 0.7 * ones(1, size(F,1));

dF                  = F - bsxfun(@times, Fneu, coefNeu(:));

% dF          = F - Fneu;

sd           = std(dF, [], 2);
sdN          = std(Fneu, [], 2);

sk(:, 1) = skewness(dF, [], 2);
sk(:, 2) = sd./sdN; 
sk(:, 3) = (max(dF, [], 2)-median(dF, 2))./sd;
sk(:, 4) = (prctile(dF, 95, 2)-median(dF, 2))./sd;

for j = 1:numel(stat)
    stat(j).dFstat           = sk(j,:);
    stat(j).skew             = sk(j,1);
    stat(j).std              = sk(j,2);
    stat(j).maxMinusMed      = sk(j,3);
    stat(j).top5pcMinusMed   = sk(j,4);
    stat(j).blockstarts      = [0 cumsum(ops.Nframes)];
    stat(j).iplane                 = ops.iplane;
    stat(j).neuropilCoefficient    = coefNeu(j); 
end
