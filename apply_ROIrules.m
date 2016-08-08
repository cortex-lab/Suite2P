function [stat, res] = apply_ROIrules(ops, stat0, res0, clustrules, flag)

if nargin<5
    flag = 0;
end
stat0    = get_regions(stat0, res0);

%%
if ops.splitROIs
    [stat, res] = get_validregions(stat0,res0, clustrules);
else
   stat = stat0;
   res = res0;
   
   for j = 1:length(stat)
      stat(j).igood = 1; 
   end
end



Ly = res0.Ly;
Lx = res0.Lx;
xs = repmat(1:Lx, Ly, 1);
ys = repmat((1:Ly)', 1, Lx);
xlx         = repmat(-ceil(Lx/2):1:ceil(Lx/2), 2*ceil(Lx/2)+1, 1);
rgrid       = sqrt(xlx.^2 + xlx'.^2);
rgridsort   = sort(rgrid(:), 'ascend');

for j = 1:length(stat)
    if numel(stat(j).lambda)==numel(stat(j).ipix)
        ipix   = stat(j).ipix;
        lambda = stat(j).lambda;
        ipix = ipix(lambda>max(lambda)/4);
        
        x0 = xs(ipix); y0 = ys(ipix);
        
        rs = ((x0 - mean(x0)).^2 + (y0 - mean(y0)).^2).^.5;
        stat(j).mrs     = mean(rs);
        stat(j).npix    = numel(ipix);
        stat(j).ipix    = ipix;
        stat(j).lambda    = lambda(lambda>max(lambda)/4);
        stat(j).mrs0    = mean(rgridsort(1:stat(j).npix));
    end
end

if flag    
    if ~exist(ops.ResultsSavePath, 'dir')
        mkdir(ops.ResultsSavePath)
    end
    
    iplane = ops.iplane;
    Nk = ops.Nk;
    
    
    save(sprintf('%s/F_%s_%s_plane%d_Nk%d.mat', ops.ResultsSavePath, ...
        ops.mouse_name, ops.date, iplane, Nk),  'ops', 'res', 'stat', 'stat0', 'res0', 'clustrules')
end
