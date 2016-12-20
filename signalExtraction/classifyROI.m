function  stat                         = classifyROI(stat, clustrules)

Mrs = [stat.mrs]./[stat.mrs0];
npix = [stat.npix];

iscell = Mrs <clustrules.Compact & ...
    npix<clustrules.MaxNpix & npix>clustrules.MinNpix;

for j = 1:length(stat)
   stat(j).iscell = iscell(j); 
end
