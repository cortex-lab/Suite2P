function add_deconvolution(ops, db)
% if neuropil was obtained by the old "model" method (pre 28.06.16)
% then set flag to 'old_model'
ops = build_ops3(db, ops);

% load the initialization of the kernel    
load(fullfile(ops.toolbox_path, 'SpikeDetection\kernel.mat'));

warning('ops0.imageRate now represents the TOTAL frame rate of the recording over all planes. This warning will be disabled in a future version. ')

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
    
    if isfield('dat', 'dat')
        dat = dat.dat; % just in case trying to load processed files
    end
    
    % overwrite fields of ops with those saved to file
    ops = addfields(ops, dat.ops);
    
    isroi = [dat.stat.mrs]./[dat.stat.mrs0]<dat.clustrules.Compact & ...
        [dat.stat.npix]>dat.clustrules.MinNpix & [dat.stat.npix]<dat.clustrules.MaxNpix;
    
    Ff = [];
    Fneu = [];
    for j = 1:numel(dat.Fcell)
        Ff   = cat(1, Ff, dat.Fcell{j}(isroi, :)');
        Fneu = cat(1,Fneu, dat.FcellNeu{j}(isroi, :)');        
    end
    
    if mean(sign(dat.FcellNeu{1}(:)))<0
        % then it means this was processed with old "model" option
        Fneu = -Fneu;
        Ff   = Ff + Fneu;
    end
    
    dcell = run_deconvolution3(ops, Ff, Fneu, kernel);
    dat.cl.isroi = isroi;    
    dat.cl.dcell = dcell;
    
    save(fpath, '-struct', 'dat')
end
%
