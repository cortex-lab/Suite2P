%%
% addpath('\\zserver\Lab\Share\Marius\TwoPhotonSuiteDev')

% make database to run in batch
% make_db_MK
make_db_5plane_red
% make_db_axons
% make_db_ss
% make_db_1plane_red
% make_db_20plane_1x;
% make_db_badreg;
% make_db_1x;

ops0.useGPU                 = 1; % if you can use a GPU in matlab this accelerate registration approx 3 times
ops0.doRegistration         = 1;
% root paths for files and temporary storage (ideally an SSD drive. my SSD is C)
ops0.RegFileTiffLocation    = []; %'D:/DATA/'; % leave empty to NOT save registered tiffs
ops0.RegFileRoot            = 'C:/DATA/tempreg';
ops0.LoadRegMean   			= 0; % 

ops0.getROIs                = 1;
ops0.getSVDcomps            = 1;
ops0.nSVD                   = 1000; % how many SVD components to keep

ops0.CopyDataLocally        = 0;
ops0.useImRead              = 0; % imread works faster from a local drive
ops0.TempStorage            = 'C:/DATA/'; % copy data locally first
ops0.ResultsSavePath        = 'D:/DATA/F';
ops0.PhaseCorrelation       = 1; % set to 0 for non-whitened cross-correlation
ops0.SubPixel               = Inf; % 2 is alignment by 0.5 pixel, Inf is the exact number from phase correlation

ops0.showTargetRegistration = 1;
ops0.RootStorage            = '//zserver4/Data/2P';
ops0.ShowCellMap            = 1;
ops0.DeleteBin              = 1; % set to 1 for batch processing on a limited hard drive
ops0.DeleteRawOnline        = 0; % set to 1 for batch processing on an (even more) limited hard drive

% these are modifiable settings for classifying ROIs post-clustering
clustrules.MaxNpix                          = 200; % important
clustrules.MinNpix                          = 20; % important
clustrules.Compact                          = 1.3; % important
clustrules.parent.minPixRelVar              = 1/10;
clustrules.parent.PixelFractionThreshold    = 0.5; % 1/20;
clustrules.parent.MaxRegions                = 10;

% the following settings shouldn't need to be adjusted
ops0.NavgFramesSVD          = 3000; % how many (pooled) frames to do the SVD based on
ops0.Nk0                    = 100;  % how many clusters to start with
ops0.Nk                     = 500;  % how many clusters to end with
ops.nSVDforROI              = 1000;
ops0.niterclustering        = 30;   % how many iterations of clustering

ops0.showTargetRegistration = 1;
ops0.NimgFirstRegistration  = 1000; 
ops0.RegPrecision           = 'int16';
ops0.RawPrecision           = 'int16';
ops0.NiterPrealign          = 10;

% ops0.sig                    = 0.25;  % spatial smoothing constant
% (encourages colocalized clusters) OBSOLETE

%%
for iexp = [2] %2:length(db)          
    % copy files from zserver
    ops = build_ops(db(iexp), ops0);
    if ops.CopyDataLocally
        copy_from_zserver(ops);
    end
    
    for iplane = 5 %ops.planesToProcess
        ops         = build_ops(db(iexp), ops0); % reset ops for each plane
        if ops.useGPU
            gpuDevice(1);   % reset GPU at each dataset
        end
        ops.iplane  = iplane;       
        ops         = reg2P(ops);  % do registration

%         fname = sprintf('%s/%s/%s/SVDmaskForROI_%s_%s_plane%d.mat', ops.ResultsSavePath, ...
%             ops.mouse_name,ops.date, ops.mouse_name, ops.date, iplane);
%         load(fname);
%         ops.Nk0 = 5000;
%         ops.Nk  = 2500;        
        if ops.getSVDcomps
            ops    = get_svdcomps(ops);
        end

        if ops.getROIs 
            [ops, U, Sv]        = get_svdForROI(ops);
            [ops, stat0, res0]  = fast_clustering(ops, reshape(U, [], size(U,3)), Sv);        
            [stat, res]         = apply_ROIrules(ops, stat0, res0, clustrules);
            Fcell               = get_signals(ops, iplane);
        end
        
        if ops.DeleteBin
            delete(ops.RegFile);        % delete temporary bin file
        end
    end
    
    % remove raw tiff directories
    if ops.CopyDataLocally && ~strfind(ops.TempStorage, 'zserver')
        % check again if this location is on zserver
        if strcmp(ops.TempStorage(1), '\') || strcmp(ops.TempStorage(1), '/')
            error('you are trying to remove a file from a network location')
        else
            rmdir(fullfile(ops.TempStorage, ops.mouse_name), 's');
        end
    end
    % clean up
    fclose all;        
end
%%
