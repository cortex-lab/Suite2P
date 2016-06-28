function add_deconvolution(ops, db, flag)
% if neuropil was obtained by the old "model" method (pre 28.06.16)
% then set flag to 'old_model'
ops = build_ops3(db, ops);

if nargin<2
    flag = 'standard';
end
% load the initialization of the kernel    
load(fullfile(ops0.toolbox_path, 'SpikeDetection\kernel.mat'));

for i = 4 %1:length(ops.planesToProcess)
    iplane  = ops.planesToProcess(i);
    
    fpath = sprintf('%s/F_%s_%s_plane%d_Nk%d.mat', ops.ResultsSavePath, ...
        ops.mouse_name, ops.date, iplane, ops.Nk);
    dat = load(fpath);
    if isfield('dat', 'dat')
        dat = dat.dat; % just in case trying to load processed files
    end
    
    isroi = [dat.stat.mrs]./[dat.stat.mrs0]<dat.clustrules.Compact & ...
        [dat.stat.npix]>dat.clustrules.MinNpix & [dat.stat.npix]<dat.clustrules.MaxNpix;
    
    Ff = [];
    Fneu = [];
    for j = 1:numel(dat.Fcell)
        Ff   = cat(1, Ff, dat.Fcell{1}(isroi, :)');
        Fneu = cat(1,Fneu, dat.FcellNeu{1}(isroi, :)');        
    end
    switch flag
        case 'old_model'
            Fneu = -Fneu;
            Ff   = Ff + Fneu;
    end
    
    dat.stat(isroi) = run_deconvolution3(ops, Ff, Fneu, kernel, dat.stat(isroi));
    dat.cl.isroi = isroi;    
    
    save('fpath', '-struct', 'dat')
end
%
