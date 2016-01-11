%% add in red channel information
addpath(genpath('\\zserver\Code\Register\'))

for iexp = 7 %2:length(db)
    if ~isempty(db(iexp).expred)
        ops = build_ops(db(iexp), ops0);
           
        for iplane = ops.planesToProcess
            
            keyboard;
            load(sprintf('%s/F_%s_%s_plane%d_Nk%d.mat', ops.ResultsSavePath,...
                ops.mouse_name, ops.date, iplane, ops.Nk))
            
            ops = getPVimage(ops);

        end
%         save(sprintf('%s/%s/%s/F_%s_%s_plane%d_Nk%d.mat', ops.ResultsSavePath, ops.mouse_name, ops.date, ...
%             ops.mouse_name, ops.date, iplane, ops.Nk), 'ops', 'res', 'stat', 'stat0', 'res0', 'Fcell')
    end
end
%% add in allinfo, pupil and running information
addpath(genpath('C:\CODE\MariusBox\Primitives'))
% iplane = 1;

for iexp = [10]
    ops.mouse_name = db(iexp).mouse_name;
    ops.date       = db(iexp).date ;
    load(sprintf('D:/DATA/F/F_%s_%s_plane%d.mat', ops.mouse_name, ops.date, iplane))
    ops.loadallinfo = 1;
    ops.loadpupil   = 1;    
    
    
    for isess = 1:numel(ops.Nframes)
        ops = collect_stimulus_stats(ops, isess);
    end
    
    save(sprintf('D:/DATA/F/F_%s_%s_plane%d.mat', ops.mouse_name, ops.date, iplane), ...
            'Fcell', 'ops', 'res', 'stat');
end

%%



