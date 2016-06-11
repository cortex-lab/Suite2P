%%
cd('D:\CODE\MariusBox\runSuite2P') % start this code in the directory with make_db
make_db_example;

toolbox_path = 'D:\CODE\GitHub\Suite2P';
if exist(toolbox_path, 'dir')
	addpath(toolbox_path) % add local path to the toolbox
else
	error('toolbox_path does not exist, please change toolbox_path');
end
ops0.clustModel  = 'neuropil'; % standard or neuropil
ops0.neuropilSub = 'surround'; % none, surround or model

ops0.useGPU                 = 0; % if you can use a GPU in matlab this accelerate registration approx 3 times
ops0.doRegistration         = 1;

% root paths for files and temporary storage (ideally an SSD drive. my SSD is C)
ops0.RootStorage            = '//zserver4.ioo.ucl.ac.uk/Data/2P';
ops0.temp_tiff              = 'C:/DATA/temp.tif'; % copy data locally first
ops0.RegFileRoot            = 'C:/DATA/';  % location for binary file
ops0.DeleteBin              = 1; % set to 1 for batch processing on a limited hard drive
ops0.ResultsSavePath        = 'D:/DATA/F';
ops0.RegFileTiffLocation    = []; %'D:/DATA/'; % leave empty to NOT save registered tiffs


ops0.showTargetRegistration = 1;
ops0.PhaseCorrelation       = 1; % set to 0 for non-whitened cross-correlation
ops0.SubPixel               = Inf; % 2 is alignment by 0.5 pixel, Inf is the exact number from phase correlation
ops0.registrationUpsample   = 1;% upsampling factor during registration, 1 for no upsampling is fastest, 2 gives better subpixel accuracy
ops0.NimgFirstRegistration  = 500; 

ops0.ShowCellMap            = 1;
ops0.getROIs                = 1;
ops0.Nk0                    = 1300;  % how many clusters to start with
ops0.Nk                     = 650;  % how many clusters to end with
ops0.sig                    = 0.5;  % spatial smoothing length for clustering; encourages localized clusters
ops0.nSVDforROI             = 1000;
ops0.nSVD                   = 1000; % how many SVD components to keep
ops0.NavgFramesSVD          = 5000; % how many (binned) timepoints to do the SVD based on
ops0.getSVDcomps            = 0; % whether to save SVD components to disk for later processing

ops0.niterclustering        = 50;   % how many iterations of clustering


%%
for iexp = 1 %:length(db)        
     run_pipeline(db(iexp), ops0, clustrules);
end
%%
