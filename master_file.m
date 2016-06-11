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

ops0.useGPU                 = 1; % if you can use a GPU in matlab this accelerate registration approx 3 times
ops0.doRegistration         = 1;
% root paths for files and temporary storage (ideally an SSD drive. my SSD is C)
ops0.RegFileTiffLocation    = []; %'D:/DATA/'; % leave empty to NOT save registered tiffs
ops0.RegFileRoot            = 'C:/DATA/tempreg';

ops0.getROIs                = 1;
ops0.getSVDcomps            = 0;
ops0.nSVD                   = 1000; % how many SVD components to keep

ops0.temp_tiff              = 'C:/DATA/temp.tiff'; % copy data locally first
ops0.ResultsSavePath        = 'D:/DATA/F';
ops0.PhaseCorrelation       = 1; % set to 0 for non-whitened cross-correlation
ops0.SubPixel               = Inf; % 2 is alignment by 0.5 pixel, Inf is the exact number from phase correlation

ops0.showTargetRegistration = 1;
ops0.RootStorage            = '//zserver4.ioo.ucl.ac.uk/Data/2P';
ops0.ShowCellMap            = 1;
ops0.DeleteBin              = 1; % set to 1 to delete bin file after processing

% these are modifiable settings for classifying ROIs post-clustering
clustrules.MaxNpix                          = 500; % important
clustrules.MinNpix                          = 30; % important
clustrules.Compact                          = 2; % important
clustrules.parent.minPixRelVar              = 1/10;
clustrules.parent.PixelFractionThreshold    = 0.5; % 1/20;
clustrules.parent.MaxRegions                = 10;

% the following settings shouldn't need to be adjusted
ops0.NavgFramesSVD          = 5000; % how many (pooled) frames to do the SVD based on
ops0.Nk0                    = 1300;  % how many clusters to start with
ops0.Nk                     = 650;  % how many clusters to end with
ops0.nSVDforROI             = 1000;
ops0.niterclustering        = 50;   % how many iterations of clustering

ops0.NimgFirstRegistration  = 500;
ops0.RegPrecision           = 'int16';
ops0.RawPrecision           = 'int16';
ops0.NiterPrealign          = 10;

ops0.LoadRegMean   			= 0; %
ops0.nimgbegend             = 250; % how many frames to average at the beginning and end of each experiment

% parameters for signal and neuropil calculation
ops1.inNeurop   = 3; % inner diameter of neuropil mask
ops1.outNeurop  = 30; % outer diameter of neuropil mask
ops1.microID    = 'b2'; % ID of microscope (b: B-scope, b2: Bergamo2, m: MOM)
ops1.useSVD     = 1; % uses all SVD components to calculate signal and
                     % neuropil instead of bin-file
ops1.getSignal  = 1; % to extract neural signal
ops1.getNeuropil= 1; % to extract neuropil
% NOTE: to save time, first run detected ROIs through the gui to select
% "good ROIs", then extract signals and neuropil only on those ROIs;
% to do so, set .getSignal=0 and .getNeuropil=0 in the master_file, then
% run get_signals_and_neuropi with the following parameters
% ops1.processed  = 1; % to load processed file
% ops1.newFile    = 1; % to create new file ending _new.mat, otherwise
%                        processed file is overwritten
% ops1.useSVD     = 1; % only choose 0 if the temporary bin-file still
%                        contains the data of the current experiment

ops0.sig                    = 0;  % spatial smoothing constant
% (encourages colocalized clusters) OBSOLETE

% parameters for signal and neuropil calculation
ops1.inNeurop   = 3; % inner diameter of neuropil mask
ops1.outNeurop  = 30; % outer diameter of neuropil mask
ops1.microID    = 'b2'; % ID of microscope (b: B-scope, b2: Bergamo2, m: MOM)
ops1.useSVD     = 1; % uses all SVD components to calculate signal and
                     % neuropil instead of bin-file
ops1.getSignal  = 1; % to extract neural signal
ops1.getNeuropil= 1; % to extract neuropil
% NOTE: to save time, first run detected ROIs through the gui to select
% "good ROIs", then extract signals and neuropil only on those ROIs;
% to do so, set .getSignal=0 and .getNeuropil=0 in the master_file, then
% run get_signals_and_neuropi with the following parameters
% ops1.processed  = 0; % to load processed file
% ops1.newFile    = 1; % to create new file ending _new.mat, otherwise
%                        processed file is overwritten
% ops1.useSVD     = 1; % has to be SVD as bin-file is only temporary
% set other parameters accordingly


%%
for iexp = 1:length(db)        %3:length(db)
    % copy files from zserver
    run_pipeline(db(iexp), ops0, clustrules);
end
%%
