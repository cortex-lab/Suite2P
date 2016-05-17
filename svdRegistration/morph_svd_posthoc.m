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
for iplane = [2:10]
    clearvars -except iplane  clustrules
    isSVD = 1;
    isGPU = 1;
    nBest = 1;

    root = 'E:\DATA';
    make_db_MP019;
%     make_db_20plane_1x;
    iExpts = [32];
%     iExpts = [17 19 22 23 24 25];
    
    [Uall,dat] = LoadMultiDay(db,iplane,iExpts, root);
    Um = Uall(:,:,1:10,:);
%     [Uall,dat,ops] = RunShiftSVDs(Uall,dat,nBest,isGPU,isSVD);
    
    ResultsSavePath = sprintf('E://DATA//F//%s//%s', ...
            db(iExpts(nBest)).mouse_name, db(iExpts(nBest)).date);
    fname = sprintf('F_%s_%s_plane%d_Nk1500.mat', db(iExpts(nBest)).mouse_name, ...
        db(iExpts(nBest)).date, iplane);
    
%     load(fullfile(ResultsSavePath,  fname))
    
%    load('E:\DATA\F\M150824_MP019\ops_all.mat')
    root = 'E:\DATA\F\M150824_MP019\multiDay';
    fname = sprintf('plane%d.mat', iplane);
    
    load(fullfile(root, fname)) 
    
%     refSVD = single(Um(:,:,:,nBest));
    for nD = [1:length(dat)]
        A = single(Um(:,:,:,nD));
        [res0,pixInv] = ShiftMasks(res,dat{nD}.ops, ops,A,refSVD,isGPU,isSVD);        
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



