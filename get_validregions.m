function [stat, res, stat0] = get_validregions(stat0, res0, clustrules)

ix = 0;

iclust = res0.iclust;
% ipixbad = true(size(res0.M));
for k = 1:length(stat0)
    
    if ~isempty(stat0(k).region)
        
%         vs = [stat0(k).region.V];
%         criterion = sum(vs > mean(vs(:)) * clustrules.parent.PixelFractionThreshold);
        if 1 %criterion<clustrules.parent.MaxRegions
            
            igood = [stat0(k).region.mrs]./[stat0(k).region.mrs0] <clustrules.Compact;
            igood = igood & [stat0(k).region.npix]>clustrules.MinNpix & ...
                [stat0(k).region.npix]<clustrules.MaxNpix;
%             igood = [stat0(k).region.npix]>clustrules.MinNpix;
            
            igood = find(igood);
            for j = 1:length(igood)
                ix = ix+1;
                stat(ix)             = stat0(k).region(igood(j));
                
                iclust(stat(ix).ipix) = ix + numel(stat0);
                
                ismmb = ismember(stat0(k).ipix, stat(ix).ipix);
                stat0(k).ipix(ismmb) = [];
            end
        end
    end
    stat0(k).mrs   = Inf;
    
end

for k = 1:length(stat0)
    stat0(k).igood = 0;
end
for k = 1:length(stat)
    stat(k).igood = 1;
end

stat(1).region = [];
stat0(1).V = [];

stat0(1).parent = [];

res = res0;
res.iclust = iclust;

% stat(1).ipix = find(res.iclust(:)==1);
% stat(1).mrs = Inf;
% stat(1).mrs0 = Inf;
% stat(1).npix = numel(stat(1).ipix);
% stat(1).lambda = res.lambda(stat(1).ipix);
% stat(1).V = Inf;

stat = [stat0 stat];
