% check out the README file first
addpath('D:\CODE\MariusBox\runSuite2P') % add the path to your make_db file
make_db_example;

toolbox_path = 'D:\CODE\GitHub\Suite2P';
if exist(toolbox_path, 'dir')
	addpath(toolbox_path) % add local path to the toolbox
else
	error('toolbox_path does not exist, please change toolbox_path');
end
ops0.clustModel             = 'neuropil'; % standard or neuropil
ops0.neuropilSub            = 'model'; % none, surround or model

ops0.useGPU                 = 0; % if you can use an Nvidia GPU in matlab this accelerate registration approx 3 times. You only need the Nvidia drivers installed (not CUDA).
ops0.doRegistration         = 1; % skip if the data is already registered

% root paths for files and temporary storage (ideally an SSD drive. my SSD is C:/)
ops0.RootStorage            = '/server/Data/2P'; % Suite2P assumes a folder structure, check out README file
ops0.temp_tiff              = 'C:/DATA/temp.tif'; % copy data locally first
ops0.RegFileRoot            = 'C:/DATA/';  % location for binary file
ops0.DeleteBin              = 0; % set to 1 for batch processing on a limited hard drive
ops0.ResultsSavePath        = 'D:/DATA/F'; % a folder structure is created inside
ops0.RegFileTiffLocation    = []; %'D:/DATA/'; % leave empty to NOT save registered tiffs

ops0.showTargetRegistration = 1;
ops0.PhaseCorrelation       = 1; % set to 0 for non-whitened cross-correlation
ops0.SubPixel               = Inf; % 2 is alignment by 0.5 pixel, Inf is the exact number from phase correlation
ops0.registrationUpsample   = 1;% upsampling factor during registration, 1 for no upsampling is much faster, 2 may give better subpixel accuracy
ops0.NimgFirstRegistration  = 500; % number of images to include in the first registration pass 

ops0.ShowCellMap            = 1;
ops0.getROIs                = 1;
ops0.Nk0                    = 1300;  % how many clusters to start with
ops0.Nk                     = 650;  % how many clusters to end with
ops0.sig                    = 0.5;  % spatial smoothing length in pixels; encourages localized clusters
ops0.nSVDforROI             = 1000;
ops0.nSVD                   = 1000; % how many SVD components to keep
ops0.NavgFramesSVD          = 5000; % how many (binned) timepoints to do the SVD based on
ops0.getSVDcomps            = 0; % whether to save SVD components to disk for later processing

ops0.niterclustering        = 50;   % how many iterations of clustering

clustrules.diameter         = 10; % expected diameter of cells (used for scaling)

db0 = db;
%%
for iexp = 1 %:length(db)        
     run_pipeline(db(iexp), ops0, clustrules);
end
%%
