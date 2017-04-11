function add_deconvolution(ops, db)
ops = build_ops3(db, ops);
ops0 = ops;


% warning('ops0.imageRate now represents the TOTAL frame rate of the recording over all planes. This warning will be disabled in a future version. ')

for i = 1:length(ops.planesToProcess)
    iplane  = ops.planesToProcess(i);
    
    fpath = sprintf('%s/F_%s_%s_plane%d_proc.mat', ops.ResultsSavePath, ...
        ops.mouse_name, ops.date, iplane);
    if exist(fpath, 'file')
        load(fpath);
    else
        fpath = sprintf('%s/F_%s_%s_plane%d.mat', ops.ResultsSavePath, ...
            ops.mouse_name, ops.date, iplane);
        dat = load(fpath);
    end
    
    if isfield(dat, 'dat')
        dat = dat.dat; % just in case...
    end
    
    % overwrite fields of ops with those saved to file
    ops = addfields(ops, dat.ops);
    
    % set up options for deconvolution
    ops.imageRate    = getOr(ops0, {'imageRate'}, 30); % total image rate (over all planes)
    ops.sensorTau    = getOr(ops0, {'sensorTau'}, 2); % approximate timescale in seconds
    ops.sameKernel   = getOr(ops0, {'sameKernel'}, 1); % 1 for same kernel per plane, 0 for individual kernels (not recommended)
    ops.sameKernel   = getOr(ops0, {'sameKernel'}, 1);
    ops.maxNeurop    = getOr(ops0, {'maxNeurop'}, Inf); % maximum neuropil coefficient (default no max)
    ops.recomputeKernel    = getOr(ops0, {'recomputeKernel'}, 0);
    
    
    fprintf('Spike deconvolution, plane %d... \n', iplane)
    
    % construct Ff and Fneu
    Ff = [];
    Fneu = [];
    for j = 1:numel(dat.Fcell)
        Ff   = cat(1, Ff, dat.Fcell{j}');
        Fneu = cat(1,Fneu, dat.FcellNeu{j}');
    end
    
    ops.fs                  = ops.imageRate/ops.nplanes;
    ops.lam                 = getOr(ops0, 'lam', 0);
    ops.estimateNeuropil    = getOr(ops0, 'estimateNeuropil', 1);
    ops.runningBaseline     = 0;

    [sp, ~, coefs,~, sd, ops] = wrapperDECONV(ops, Ff, Fneu);
    
    stat = dat.stat;
    for j = 1:size(Ff,2)
        stat(j).neuropilCoefficient = coefs(j);
        stat(j).noiseLevel          = sd(j);  
    end
    
    nCum = 0;
    for j = 1:length(dat.Fcell)
       dat.sp{j} = sp(nCum + [1:size(dat.Fcell{j},2)], :)'; 
       nCum = nCum + size(dat.Fcell{j},2);
    end
    
    
    % compute neuropil coefficient and run deconvolution on
    % neuropil-corrected trace
    %     stat = run_deconvolution3(ops, dat);
    
    
    dat.stat = stat;
    
    save(fpath, '-struct', 'dat')
    
end
%
