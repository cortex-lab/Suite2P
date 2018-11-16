function stat = getNonOverlapROIs(stat, Ly, Lx)


Mask = zeros(Ly, Lx);

for k = 1:numel(stat)
   Mask(stat(k).ipix) = Mask(stat(k).ipix) + 1;
end

for k = 1:numel(stat)
   stat(k).isoverlap = Mask(stat(k).ipix)>1; 
end
