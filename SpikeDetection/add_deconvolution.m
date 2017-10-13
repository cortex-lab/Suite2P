% adds deconvolved traces (sp) to files
% by default it optimizes the neuropil coefficients
% max neuropil coefficients set in ops.maxNeurop

function add_deconvolution(ops, db, iplane0)
ops = build_ops3(db, ops);
ops0 = ops;

% ops.deconvType = 'OASIS' (default) or 'L0'. See
% http://www.biorxiv.org/content/early/2017/06/27/156786 for more info. 

% warning: ops0.imageRate now represents the TOTAL frame rate of the recording over all planes. This warning will be disabled in a future version.

for i =  1:length(ops.planesToProcess)    
    iplane  = ops0.planesToProcess(i);
    
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
    
    % if z-drift computed, apply correction
    for j = 1:length(dat.Fcell)
        Fcell{j}    = dat.Fcell{j};
        FcellNeu{j} = dat.FcellNeu{j};
    end
    
    % overwrite fields of ops with those saved to file
    if getOr(ops, 'newops', 0)
        ops = addfields(dat.ops, ops0); % overwrite saved ops with new
    else
        ops = addfields(ops0, dat.ops);% overwrite new ops with saved
    end
    
    % set up options for deconvolution
    ops.imageRate    = getOr(ops0, {'imageRate'}, 30); % total image rate (over all planes)
    ops.sensorTau    = getOr(ops0, {'sensorTau'}, 2); % approximate timescale in seconds
    ops.sameKernel   = getOr(ops0, {'sameKernel'}, 1); % 1 for same kernel per plane, 0 for individual kernels (not recommended)
    ops.sameKernel   = getOr(ops0, {'sameKernel'}, 1);
    ops.maxNeurop    = getOr(ops0, {'maxNeurop'}, Inf); % maximum neuropil coefficient (default no max)
    ops.recomputeKernel    = getOr(ops0, {'recomputeKernel'}, 0);    
    ops.deconvNeuropil = getOr(ops0, {'deconvNeuropil'}, 0);   % whether to deconvolve the neuropil as well
    
    fprintf('Spike deconvolution, plane %d... \n', iplane)
    
    % construct Ff and Fneu
    Ff = [];
    Fneu = [];
    for j = 1:numel(Fcell)
        Ff   = cat(1, Ff, Fcell{j}');
        Fneu = cat(1,Fneu, FcellNeu{j}');
    end
    
    ops.fs                  = ops.imageRate/ops.nplanes;
    ops.estimateNeuropil    = getOr(ops0, 'estimateNeuropil', 1);

    if ops.estimateNeuropil
        [sp, ~, coefs,~, sd, ops, baselines] = wrapperDECONV(ops, Ff, Fneu);
    else
        if isfield(dat.stat, 'neuropilCoefficient')
            coefs = [dat.stat.neuropilCoefficient];
        else
            coefs = .8 * ones(size(Ff,2), 1);
        end
    end
    
    if ops.deconvNeuropil
        ops.estimateNeuropil = 0;
        spNeu = wrapperDECONV(ops, Fneu);
    end
    
    stat = dat.stat;
    for j = 1:size(Ff,2)
        stat(j).neuropilCoefficient = coefs(j);
    end
    
    ops.SpkDcnvAlignToPlane = getOr(ops, 'SpkDcnvAlignToPlane', 0);
    if getOr(ops, 'runningBaseline', 0)
        Ff = Ff - bsxfun(@times, Fneu, [stat.neuropilCoefficient]);
        if ops.SpkDcnvAlignToPlane
            ds = (ops.SpkDcnvAlignToPlane - iplane)/ops.nplanes;
            Ff    = register_F(Ff, ds);
        end
        [sp, ~, ~,~, sd, ops, baselines] = wrapperDECONV(ops, Ff);
    end
    
    for j = 1:size(Ff,2)
        stat(j).noiseLevel          = sd(j);  
        stat(j).baseline           = baselines(j);
    end
    
    nCum = 0;
    for j = 1:length(dat.Fcell)
       dat.sp{j} = sp(nCum + [1:size(dat.Fcell{j},2)], :)'; 
       if ops.deconvNeuropil
           dat.spNeu{j} = spNeu(nCum + [1:size(dat.Fcell{j},2)], :)';
       end
       nCum = nCum + size(dat.Fcell{j},2);
    end
    
    % compute neuropil coefficient and run deconvolution on
    % neuropil-corrected trace
    %     stat = run_deconvolution3(ops, dat);
    
    dat.stat = stat;    
        
    save(fpath, '-struct', 'dat')    
end
%
