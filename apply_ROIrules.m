function [stat, res] = apply_ROIrules(ops, stat0, res0, clustrules)

stat0    = get_regions(stat0, res0, clustrules);

%%
if 1
    [stat, res] = get_validregions(stat0,res0, clustrules);
else
   stat = stat0;
   res = res0;
end


if ~exist(ops.ResultsSavePath, 'dir')
    mkdir(ops.ResultsSavePath)
end

iplane = ops.iplane;
Nk = ops.Nk;


save(sprintf('%s/F_%s_%s_plane%d_Nk%d.mat', ops.ResultsSavePath, ...
    ops.mouse_name, ops.date, iplane, Nk),  'ops', 'res', 'stat', 'stat0', 'res0')

