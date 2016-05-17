% master_fileMP
addpath('C:\CODE\GitHub\Suite2P\svdRegistration')
addpath('C:\CODE\GitHub\Suite2P\')
addpath('D:\CODE\MariusBox\runSuite2P')
addpath('C:\CODE\GitHub\Suite2P\multiDay')
%
% these are modifiable settings for classifying ROIs post-clustering
clustrules.MaxNpix                          = 200; % important
clustrules.MinNpix                          = 10; % important
clustrules.Compact                          = 2; % important
clustrules.parent.minPixRelVar              = 0;
clustrules.parent.PixelFractionThreshold    = 0; % 1/20;
clustrules.parent.MaxRegions                = Inf;

%

for iplane = 5 %2:10
    clearvars -except iplane  clustrules
    isSVD = 1;
    isGPU = 1;
    nBest = 5;

    root = 'E:\DATA';
    make_db_20plane_1x;
    iExpts = [17 19 22 23 24 25];
    
    [Uall,dat] = LoadMultiDay(db,iplane,iExpts, root);
    Um = Uall(:,:,1:10,:);
    
    % shift svds to mean image of nBest and multiply Uall by Sv
    [Uall,dat,ops] = RunShiftSVDs(Uall,dat,nBest,isGPU, isSVD);
    
    %
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
    B = single(Um(:,:,:,nBest));
    for nD = [1:length(dat)]
        A = single(Um(:,:,:,nD));
        [res0,pixInv] = ShiftMasks(res,ops,A,B,isGPU,isSVD);        
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
    for nD = 1:numel(dat)
        opsi = dat{nD}.ops;
        opsi.yrange = dat{nD}.res.yrange;
        opsi.xrange = dat{nD}.res.xrange;
        
        opsi.RegFile = sprintf('E://DATA//%s//%s//plane%d.bin', ...
            db(iExpts(nD)).mouse_name, db(iExpts(nD)).date, iplane);
        opsi.ResultsSavePath = sprintf('E://DATA//F//%s//%s', ...
            db(iExpts(nD)).mouse_name, db(iExpts(nD)).date);
        
        resT = dat{nD}.res;
        resT.iclust(resT.iclust(:)==0) = 1;
        statT = dat{nD}.stat;
        get_signals_NEUmodel_multi(statT, resT, opsi, iplane);
    end
end



