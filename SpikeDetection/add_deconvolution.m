function add_deconvolution(ops, db, clustrules)
ops = build_ops3(db, ops);
ops0 = ops;

clustrules = get_clustrules(clustrules);

try
    ppool = gcp ;
catch
end

% warning('ops0.imageRate now represents the TOTAL frame rate of the recording over all planes. This warning will be disabled in a future version. ')

for i = 1:length(ops.planesToProcess)
    iplane  = ops.planesToProcess(i);
    
    fpath = sprintf('%s/F_%s_%s_plane%d_Nk%d_proc.mat', ops.ResultsSavePath, ...
        ops.mouse_name, ops.date, iplane, ops.Nk);
    if exist(fpath, 'file')
        load(fpath);
    else
        fpath = sprintf('%s/F_%s_%s_plane%d_Nk%d.mat', ops.ResultsSavePath, ...
            ops.mouse_name, ops.date, iplane, ops.Nk);
        dat = load(fpath);
    end
    
    if isfield(dat, 'dat')
        dat = dat.dat; % just in case trying to load processed files
    end
    
    
    % overwrite fields of ops with those saved to file
    ops = addfields(ops, dat.ops);
    if ~isfield(dat, 'clustrules');
        dat.clustrules = clustrules;
    end
    
    % set up options for deconvolution
    ops.imageRate    = getOr(ops0, {'imageRate'}, 30); % total image rate (over all planes)
    ops.sensorTau    = getOr(ops0, {'sensorTau'}, 2); % approximate timescale in seconds
    ops.sameKernel   = getOr(ops0, {'sameKernel'}, 1); % 1 for same kernel per plane, 0 for individual kernels (not recommended)
    ops.sameKernel   = getOr(ops0, {'sameKernel'}, 1);
    ops.maxNeurop    = getOr(ops0, {'maxNeurop'}, Inf);
    ops.recomputeKernel    = getOr(ops0, {'recomputeKernel'}, 0);
    
    
    fprintf('Spike deconvolution, plane %d... \n', iplane)
    % split data into batches
    [dcell, isroi] = run_deconvolution3(ops, dat);
    
    dat.cl.isroi = isroi;
    dat.cl.dcell = dcell;
    
    save(fpath, '-struct', 'dat')
end
%
