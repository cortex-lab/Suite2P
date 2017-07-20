function  run_pipeline(db, ops0)

% ops0.TileFactor (or db(iexp).TileFactor) can be set to multiply the number of default tiles for the neuropil

ops0.nimgbegend                     = getOr(ops0, {'nimgbegend'}, 0);
ops0.splitROIs                      = getOr(ops0, {'splitROIs'}, 1);
ops0.LoadRegMean                    = getOr(ops0, {'LoadRegMean'}, 0);
ops0.NiterPrealign                  = getOr(ops0, {'NiterPrealign'}, 10);
ops0.registrationUpsample           = getOr(ops0, {'registrationUpsample'}, 1);  % upsampling factor during registration, 1 for no upsampling is much faster, 2 may give better subpixel accuracy
ops0.getROIs                        = getOr(ops0, {'getROIs'}, 1);   % whether to run the optimization
ops0.getSVDcomps                    = getOr(ops0, {'getSVDcomps'}, 0);   % whether to save SVD components to disk for later processing
ops0.nSVD                           = getOr(ops0, {'nSVD'}, 1000);   % how many SVD components to save to disk
ops0.signalExtraction               = getOr(ops0, 'signalExtraction', 'raw');
ops0.interpolateAcrossPlanes        = getOr(ops0, 'interpolateAcrossPlanes', 0);
ops0.maxNeurop                      = getOr(ops0, 'maxNeurop', 1.5);

% neuropil options for surround neuropil computation
if strcmp(ops.signalExtraction, 'surround')
    ops.innerNeuropil = getOr(ops, 'innerNeuropil', 1);
    ops.outerNeuropil = getOr(ops, 'outerNeuropil', Inf);
    if isinf(ops.outerNeuropil)
        ops.minNeuropilPixels = getOr(ops, 'minNeuropilPixels', 400);
        ops.ratioNeuropil     = getOr(ops, 'ratioNeuropil', 5);
    end
end


ops                                 = build_ops3(db, ops0);
if isfield(ops, 'numBlocks') && ~isempty(ops.numBlocks)
    if numel(ops.numBlocks) == 1
        ops.numBlocks = [ops.numBlocks 1];
    end
    if sum(ops.numBlocks) > 2
        ops.nonrigid               = 1;
    end
end
ops.nonrigid                       = getOr(ops, 'nonrigid', 0);   
ops.kriging                        = getOr(ops, 'kriging', 1);  

if ~isfield(ops, 'diameter') || isempty(ops.diameter)
    warning('you have not specified mean diameter of your ROIs')
    warning('for best performance, please set db(iexp).diameter for each experiment')
end
ops.diameter                        = getOr(ops, 'diameter', 10);
ops.clustrules.diameter             = ops.diameter;
ops.clustrules                      = get_clustrules(ops.clustrules);
%
% this loads ops1 and checks if processed binary files exist
opath = sprintf('%s/regops_%s_%s.mat', ops.ResultsSavePath, ops.mouse_name, ops.date);
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
%
% do registration if the processed binaries do not exist
if processed==0
    if ops.nonrigid
        ops1 = blockReg2P(ops);  % do non-rigid registration
    else
        ops1 = reg2P(ops);  % do registration
    end
else
    disp('already registered binary found');
end
% keyboard;
%%
for i = 1:numel(ops1) %[1:2 4:numel(ops1)]
    ops         = ops1{i};    
  
    ops.iplane  = i;
    
%     ops.ThScaling = 0.5;
    
    if numel(ops.yrange)<10 || numel(ops.xrange)<10
        warning('valid range after registration very small, continuing to next plane')
        continue;
    end
    
    if getOr(ops, {'getSVDcomps'}, 0)
        % extract and write to disk SVD comps (raw data)
        ops    = get_svdcomps(ops);
    end
    
    if ops.getROIs || getOr(ops, {'writeSVDroi'}, 0)
        % extract and/or write to disk SVD comps (normalized data)
        [ops, U, model]    = get_svdForROI(ops);
    end
        
    if ops.getROIs
        % get sources in stat, and clustering images in res
        [ops, stat, model]           = sourcery(ops,U, model);
        
        % extract dF
        switch ops.signalExtraction
            case 'raw'
                [ops, stat, Fcell, FcellNeu] = extractSignalsNoOverlaps(ops, model, stat);
            case 'regression'
                [ops, stat, Fcell, FcellNeu] = extractSignals(ops, model, stat);
            case 'surround'
                [ops, stat, Fcell, FcellNeu] = extractSignalsSurroundNeuropil(ops, stat);
        end

        % apply user-specific clustrules to infer stat.iscell
        stat                         = classifyROI(stat, ops.clustrules);
        
        save(sprintf('%s/F_%s_%s_plane%d.mat', ops.ResultsSavePath, ...
            ops.mouse_name, ops.date, ops.iplane),  'ops',  'stat',...
            'Fcell', 'FcellNeu', '-v7.3')
    end

    
    if ops.DeleteBin
        fclose('all');
        delete(ops.RegFile);        % delete temporary bin file
    end
end

% clean up
fclose all;