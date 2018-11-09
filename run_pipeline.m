function  ops1 = run_pipeline(db, ops0)

% ops0.TileFactor (or db(iexp).TileFactor) can be set to multiply the number of default tiles for the neuropil

ops0.splitROIs                      = getOr(ops0, {'splitROIs'}, 1);
ops0.LoadRegMean                    = getOr(ops0, {'LoadRegMean'}, 0);
ops0.getROIs                        = getOr(ops0, {'getROIs'}, 1);   % whether to run the optimization
ops0.getSVDcomps                    = getOr(ops0, {'getSVDcomps'}, 0);   % whether to save SVD components to disk for later processing

ops0                                = build_ops3(db, ops0);

if ~isfield(ops0, 'diameter') || isempty(ops0.diameter)
    error('you must specify approximate mean diameter of your ROIs in ops.diameter or db.diameter')    
end
if ~isfield(ops0, 'tau') || isempty(ops0.tau)
    error('you must specify approximate timescale in seconds in ops.tau or db.tau')    
end
if ~isfield(ops0, 'fs') || isempty(ops0.fs)
    error('you must specify sampling rate / plane in ops.fs or db.fs')    
end

ops0.diameter                        = getOr(ops0, 'diameter', 8*ops0.zoom);
ops0.clustrules.diameter             = ops0.diameter;
ops0.clustrules                      = get_clustrules(ops0.clustrules);

% this loads ops1 and checks if processed binary files exist
opath = sprintf('%s/regops_%s_%s.mat', ops0.ResultsSavePath, ops0.mouse_name, ops0.date);
processed = 1;
if exist(opath, 'file')
    load(opath);
    for j = 1:numel(ops1)       
       if ~exist(ops1{j}.RegFile, 'file') % check if the registered binary file exists
          processed = 0; 
       end
    end
else
    processed = 0;
end

% run reg2P if the binaries do not exist
%%%% if tiffs have already been registered, set ops.doRegistration = 0
%%%% and reg2P will just create binary file
% ops1 are the settings and values from registration
if processed==0
    ops1 = reg2P(ops0);  % do registration
else
    disp('already registered binary found');
end


%%
for i = [1:numel(ops1)]
    ops         = ops1{i};    

    % check if settings are different between ops and ops0
    % ops0 settings are chosen over ops settings
    ops         = opsChanges(ops, ops0);    
    
    ops.iplane  = i;
    
    if numel(ops.yrange)<10 || numel(ops.xrange)<10
        warning('valid range after registration very small, continuing to next plane')
        continue;
    end
            
    if ops.getROIs
        % get sources in stat, and clustering images in res
        [ops, stat]           = sourcery_inc(ops);
        
        [ops, stat, Fcell, FcellNeu] = extractF(ops, stat);
        
        % apply user-specific clustrules to infer stat.iscell
        stat                         = classifyROI(stat, ops.clustrules);
        
        save(sprintf('%s/F_%s_%s_plane%d.mat', ops.ResultsSavePath, ...
            ops.mouse_name, ops.date, ops.iplane),  'ops',  'stat',...
            'Fcell', 'FcellNeu', '-v7.3')
    end

    ops1{i} = ops;
    
    fclose('all');
    if ops.DeleteBin       
         delete(ops.RegFile);        % delete temporary bin file
    end
end

% clean up
fclose all;
