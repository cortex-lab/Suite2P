% master_fileMP
for iplane = [2:10]
    clearvars -except iplane  clustrules
    root = 'E:\DATA';
    make_db_20plane_1x;
    iExpts = [17:22];
    
    [Uall,dat] = LoadMultiDay(db,iplane,iExpts, root);
    %
    addpath('C:\CODE\GitHub\Suite2P\svdRegistration')
    
    isGPU = 1;
    nBest = 4;
    
    % shift svds to mean image of nBest and multiply Uall by Sv
    [Uall,dat,ops] = RunShiftSVDs(Uall,dat,nBest,isGPU);
    Uall = Uall(ops.yrange,ops.xrange,:,:);
    % take top 1000 components of Uall
    Uall = reshape(Uall,[],size(Uall,3)*size(Uall,4));
    [Un, Sv, ~] = svdecon(Uall'*Uall);
    Sv = diag(Sv);
    
    U = Uall * Un(:, 1:1000);
    clear Uall;
    U = normc(U);
    Sv = Sv/(numel(dat)/2);
    
    U = reshape(U, numel(ops.yrange), numel(ops.xrange), []);
    %
    [ops,stat,res] = fast_clustering_with_neuropil_multi(ops,U,Sv);
    
    %
    ops.iplane = iplane;
    [stat, res] = apply_ROIrules_multi(ops, stat, res, clustrules);
    
    
    %
    B = single(dat{nBest}.ops.mimg1);
    for nD = 1:numel(dat)
        A = single(dat{nD}.ops.mimg1);
        [res0,pixInv] = ShiftMasks(res,ops,A,B,isGPU);
        
        stat0    = get_stat(res0);
        stat1    = stat;
        
        for i = 1:length(stat0)
            stat1(i).ipix = stat0(i).ipix;
        end
        
        dat{nD}.res    = res0;
        dat{nD}.stat    = stat1;
        dat{nD}.pixInv = pixInv;
    end
    
    
    
    %
    for nD = 1:6
        ops = dat{nD}.ops;
        ops.yrange = dat{nD}.res.yrange;
        ops.xrange = dat{nD}.res.xrange;
        
        ops.RegFile = sprintf('E://DATA//%s//%s//plane%d.bin', db(iExpts(nD)).mouse_name, db(iExpts(nD)).date, iplane);
        ops.ResultsSavePath = sprintf('E://DATA//F//%s//%s', db(iExpts(nD)).mouse_name, db(iExpts(nD)).date);
        
        res = dat{nD}.res;
        res.iclust(res.iclust(:)==0) = 1;
        stat = dat{nD}.stat;
        get_signals_NEUmodel_multi(stat, res, ops, iplane);
    end
end



