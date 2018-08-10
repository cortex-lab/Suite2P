%% SET ALL DEFAULT OPTIONS HERE

% UPDATE fall 2017: non-rigid and rigid registration scripts merged; red
% channel mean image can be computed while registering green channel; red
% channel binary can be computed during green channel registration
% (ly x lx x time like green channel)

% UPDATE end-of-summer 2017: default neuropil extraction is now "surround"
% and it's very fast. Cell extraction is on the raw data (no pixel-scaling or smoothing).

% UPDATE summer 2017: default spike deconvolution changed to a customized version of
% OASIS (due to our results in this paper http://www.biorxiv.org/content/early/2017/06/27/156786). Please
% Please download the OASIS code from https://github.com/zhoupc/OASIS_matlab, and
% add the folder, with its subfolders, to your Matlab path.

% UPDATE Christmas 2016: number of clusters determined automatically, but
% do specify the "diameter" of an average cell for best results. You can do this with either
% db(iexp).diameter, or ops0.diameter

% check out the README file for detailed instructions
% **** and for more options available ****
addpath('D:\CODE\MariusBox\runSuite2P') % add the path to your make_db file

ops0.toolbox_path = 'C:\CODE\GitHub\Suite2P';
if exist(ops0.toolbox_path, 'dir')
    addpath(genpath(ops0.toolbox_path)) % add local path to the toolbox
else
    error('toolbox_path does not exist, please change toolbox_path');
end

% mex -largeArrayDims SpikeDetection/deconvL0.c (or .cpp) % MAKE SURE YOU COMPILE THIS FIRST FOR DECONVOLUTION

% overwrite any of these default options in your make_db file for individual experiments
make_db_example; % RUN YOUR OWN MAKE_DB SCRIPT TO RUN HERE

ops0.useGPU                 = gpuDeviceCount > 0; % if you can use an Nvidia GPU in matlab this accelerates registration approx 3 times. You only need the Nvidia drivers installed (not CUDA).
ops0.fig                    = 1; % turn off figure generation with 0
% ops0.diameter               = 12; % most important parameter. Set here, or individually per experiment in make_db file

% ---- root paths for files and temporary storage (ideally an SSD drive. my SSD is C:/)
% are defined through the db.
ops0.DeleteBin              = 1; % set to 1 for batch processing on a limited hard drive
% if you want to save red channel tiffs, also set ops0.REDbinary = 1

% ---- registration options ------------------------------------- %
ops0.doRegistration         = 1; % skip (0) if data is already registered
ops0.showTargetRegistration = 1; % shows the image targets for all planes to be registered
ops0.PhaseCorrelation       = 1; % set to 0 for non-whitened cross-correlation
ops0.SubPixel               = Inf; % 2 is alignment by 0.5 pixel, Inf is the exact number from phase correlation
ops0.NimgFirstRegistration  = 500; % number of images to include in the first registration pass
ops0.nimgbegend             = 0; % frames to average at beginning and end of blocks
ops0.dobidi                 = 1; % infer and apply bidirectional phase offset

% ---- cell detection options ------------------------------------------%
ops0.ShowCellMap            = 1; % during optimization, show a figure of the clusters
ops0.sig                    = 0.5;  % spatial smoothing length in pixels; encourages localized clusters
ops0.nSVDforROI             = 1000; % how many SVD components for cell clustering
ops0.NavgFramesSVD          = 5000; % how many (binned) timepoints to do the SVD based on
ops0.signalExtraction       = 'surround'; % how to extract ROI and neuropil signals:
%  'raw' (no cell overlaps), 'regression' (allows cell overlaps),
%  'surround' (no cell overlaps, surround neuropil model)

% ----- neuropil options (if 'surround' option) ------------------- %
% all are in measurements of pixels
ops0.innerNeuropil  = 1; % padding around cell to exclude from neuropil
ops0.outerNeuropil  = Inf; % radius of neuropil surround
% if infinity, then neuropil surround radius is a function of cell size
if isinf(ops0.outerNeuropil)
    ops0.minNeuropilPixels = 400; % minimum number of pixels in neuropil surround
    ops0.ratioNeuropil     = 5; % ratio btw neuropil radius and cell radius
    % radius of surround neuropil = ops0.ratioNeuropil * (radius of cell)
end

% ----- spike deconvolution and neuropil subtraction options ----- %
ops0.imageRate              = 30;   % imaging rate (cumulative over planes!). Approximate, for initialization of deconvolution kernel.
ops0.sensorTau              = 2; % decay half-life (or timescale). Approximate, for initialization of deconvolution kernel.
ops0.maxNeurop              = 1; % for the neuropil contamination to be less than this (sometimes good, i.e. for interneurons)

% ----- if you have a RED channel ---------------------- ------------%
ops0.AlignToRedChannel      = 0; % compute registration offsets using red channel
ops0.REDbinary              = 0; % make a binary file of registered red frames
% if db.expred, then compute mean red image for green experiments with red
% channel available while doing registration
ops0.redMeanImg             = 0;
% for red cell detection (identify_redcells_sourcery.m)
% redratio = red pixels inside / red pixels outside
% redcell = redratio > mean(redratio) + redthres*std(redratio)
% notred = redratio < mean(redratio) + redmax*std(redratio)
ops0.redthres               = 1.5; % the higher the thres the less red cells
ops0.redmax                 = 1; % the higher the max the more NON-red cells

db0 = db;

sm = CortexLabStorageManager(ops0);
for i = 1:length(db)
    sm.addEntry(db(i));
end

%% RUN THE PIPELINE HERE

for i = 1 %[1:sm.getNumEntries()]
    sm.selectEntry(i);
    run_pipeline(sm);

    % deconvolved data into st, and neuropil subtraction coef in stat
    add_deconvolution(sm);

    % add red channel information (if it exists)
    expred = sm.getOr('expred', []);
    if ~isempty(expred)
        % % creates mean red channel image aligned to green channel
        % % use this if you didn't get red channel during registration
        % % OR you have a separate experiment with red and green just for this
        % red_expts = ismember(db.expts, getOr(db, 'expred', []));
        % if ~ops0.redMeanImg || sum(red_expts)==0
        %     run_REDaddon_sourcery(db, ops0);
        % end

        % % identify red cells in mean red channel image
        % % fills dat.stat.redcell, dat.stat.notred, dat.stat.redprob
        % identify_redcells_sourcery(db, ops0);
    end

end
%% STRUCTURE OF RESULTS FILE
% cell traces are in dat.Fcell
% neuropil traces are in dat.FcellNeu
% manual, GUI overwritten "iscell" labels are in dat.cl.iscell
%
% stat(icell) contains all other information:
% iscell: automated label, based on anatomy
% neuropilCoefficient: neuropil subtraction coefficient, based on maximizing the skewness of the corrected trace (ICA)
% st: are the deconvolved spike times (in frames)
% c:  are the deconvolved amplitudes
% kernel: is the estimated kernel
