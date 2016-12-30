UPDATE (Christmas 2016): number of clusters determined automatically, but do specify the "diameter" of an average cell for best results. You can do this with
either db(iexp).diameter or ops.diameter. 

# Suite2p: fast, accurate and complete two-photon pipeline#

Registration, cell detection, spike extraction and manual GUI. 

Details in [http://biorxiv.org/content/early/2016/06/30/061507](http://biorxiv.org/content/early/2016/06/30/061507).

[![IMG](https://img.youtube.com/vi/xr-flH2Ow2Y/0.jpg)](https://www.youtube.com/watch?v=xr-flH2Ow2Y)

This code was written by Marius Pachitariu and members of the lab of Kenneth Harris and Matteo Carandini. It is provided here with no warranty. For support, please open an issue directly on github. 

### I. Introduction ###

This is a complete automated pipeline for processing two-photon Calcium imaging recordings. It is very simple, very fast and yields a large set of active ROIs. A GUI further provides point-and-click capabilities for refining the results in minutes. The pipeline includes the following steps


1) X-Y subpixel registration --- using a version of the phase correlation algorithm and subpixel translation in the FFT domain. If a GPU is available, this completes in 20 minutes per 1h of recordings at 30Hz and 512x512 resolution.

2) SVD decomposition --- this provides the basis for a number of pixel-level visualization tools. 

3) Cell detection --- using clustering methods in a low-dimensional space of the fluorescence activity. The pixel clustering algorithm directly provides a positive mask for each ROI identified and stands in contrast to more "black-box" methods that have been previously proposed. The mask is used to weigh the pixels of the ROI before taking the average over pixel ROIs. 

4) Manual curation --- the output of the cell detection algorithm can be visualized and further refined using the included GUI. The GUI is designed to make cell sorting a fun and enjoyable experience. 

5) Spike deconvolution --- cell and neuropil traces are further processed to obtain an estimate of spike times and spike "amplitudes". The amplitudes here are understood as proportional to the number of spikes in a burst/bin. Under low SNR conditions, the deconvolution is still useful for temporally-localizing responses, and for providing an estimate of the neuropil contamination coefficient (estimated together, see paper for how this works).

### II. Getting started ###

The toolbox runs in Matlab (+ one mex file) and currently only supports tiff file inputs. To begin using the toolbox, you will need to make local copies (in a separate folder) of two included files: master_file and make_db. It is important that you make local copies of these files, otherwise updating the repository will overwrite them (and you can lose your files). The make_db file assembles a database of experiments that you would like to be processed in batch. It also adds per-session specific information that the algorithm requires such as the number of imaged planes and channels. The master_file sets general processing options that are applied to all sessions included in make_db, UNLESS the option is over-ridden in the make_db file.  The global and session-specific options are described in detail below. 

For spike deconvolution, you need to run mex -largeArrayDims SpikeDetection/deconvL0.c (or .cpp under Linux/Mac)

### III. Input-output file paths ###

RootStorage --- the root location where the raw tiff files are  stored.

RegFileRoot --- location on local disk where to keep the registered movies in binary format. This will be loaded several times so it should ideally be an SSD drive. 

ResultsSavePath --- where to save the final results. 

DeleteBin --- deletes the binary file created to store the registered movies

RegFileTiffLocation --- where to save registered tiffs (if empty, does not save)

All of these filepaths are completed with separate subfolders per animal and experiment, specified in the make_db file. Your data MUST be stored under a file tree of the form

\RootStorage\mouse_name\session\block\*.tif(f)

The output is a struct called dat which is saved into a mat file in ResultsSavePath using the same subfolder structure, under a name formatted like F_M150329_MP009_2015-04-29_plane1_Nk650. It contains all the information collected throughout the processing, and contains the fluorescence traces in dat.F.Fcell and whether a given ROI is a cell or not in dat.F.iscell. dat.stat contains information about each ROI and can be used to recover the corresponding pixels for each ROI in dat.stat.ipix. The centroid of the ROI is specified in dat.stat as well. Here is a summary of where the important results are:

cell traces are in dat.F.Fcell
neuropil traces are in dat.F.FcellNeu
neuropil subtraction coefficient is dat.cl.dcell{i}.B(3)
baseline is dat.cl.dcell{i}.B(2)
anatomical cell criterion is in dat.cl.isroi
manual overwritten cell labels are in dat.cl.iscell
dat.cl.dcell{i}.st are the deconvolved spike times (in frames)
dat.cl.dcell{i}.c  are the deconvolved amplitudes
dat.cl.dcell{i}.kernel is the estimated kernel


### IV. Options for registration ###

showTargetRegistration --- whether to show an image of the target frame immediately after it is computed. 

PhaseCorrelation --- whether to use phase correlation (the alternative is normal cross-correlation).

SubPixel --- accuracy level of subpixel registration (10 = 0.1 pixel accuracy)

kriging --- compute shifts using kernel regression with a gaussian kernel of width 1 onto a grid of 1/SubPixel

NimgFirstRegistration --- number of randomly sampled images to do the target computation from

NiterPrealign --- number of iterations for the target computation (iterative re-alignment of subset of frames)

smooth_time_space --- convolves raw movie with a Gaussian of specified size in specified dimensions;
                      [t]: convolve in time with gauss. of std t, [t s]: convolve in time and space,
                      [t x y]: convolve in time, and in space with an ellipse rather than circle
                      
++ Block Registration (for high zoom/npixels - assumes scanning is in Y direction) ++

nonrigid --- set to 1 for non-rigid registration (or set numBlocks > 1)

numBlocks --- number of blocks dividing y-dimension of image (default is 6)

blockFrac --- percent of image to use per block (default is 1/(numBlocks-1))
             
quadBlocks --- interpolate block shifts to single line shifts (6 blocks -> 512 lines) by fitting a quadratic function (default is 1)

smoothBlocks --- if quadBlocks = 0, then smoothBlocks is the standard deviation of the gaussian smoothing kernel

++ Bidirectional scanning issues (frilly cells) taken care of automatically ++

### V. Options for cell detection ###

clustModel --- what clustering model to use; "neuropil" is described in the paper, adds a neuropil contribution to each pixel individually. "standard" is plain clustering. 

neuropilSub --- for neuropil subtraction. "none" does not output a neuropil trace; "surround" estimates the signal in an annulus around the cell after ignoring all ROIs detected during optimization (even non-cells); "model" is described in the paper, uses spatial basis functions to estimate a smoothly-varying neuropil signal over the FOV and regresses this out simultaneously with the estimation of ROI responses. 

sig --- spatial smoothing constant: smooths out the SVDs spatially. Makes ROIs more cell-shaped. 

Nk0 --- starting the algorithm with this many clusters

Nk --- final annealed number of clusters (warning, due to the structure of the annealing process, one should have Nk0<3*Nk or else the code crashes. No need to make it larger anyway). 

nSVDforROI --- how many SVD components to keep for clustering. Usually ~ the number of expected cells (Nk). 

ShowCellMap --- whether to show the clustering results as an image every 10 iterations of the clustering

niterclustering --- how many iterations of clustering

getROIs --- whether to run the ROI detection algorithm after registration

### Options for SVD decomposition ###

NavgFramesSVD --- for SVD, data has to be temporally binned. This number specifies the final number of points to be obtained after binning. In other words, datasets with many timepoints are binned in higher windows while small dataset are binned less. 

getSVDcomps --- whether to obtain and save to disk SVD components of the registered movies. Useful for pixel-level analysis and for checking the quality of the registration (residual motion will show up as SVD components). This is a separate SVD decomposition from that done for cell clustering (does not remove a running baseline of each pixel). 

nSVD --- how many SVD components to keep.

### VII. Options for spike deconvolution ###

imageRate --- imaging rate per plane. Approximate, for initialization of deconvolution kernel.  

sensorTau --- decay half-life (or timescale). Approximate, for initialization of deconvolution kernel.

maxNeurop --- neuropil contamination coef has to be less than this (sometimes good, i.e. for interneurons)

### VIII. Rules for post-clustering ROI classification ###

These options serve to compute candidate cell clusters, that can then be refined in the GUI. Clusters computed in the algorithm are split into connected regions and then classified as cell/non-cell based on the following

diameter --- % expected diameter of cells (used for 0.25 * pi/4*diam^2 < npixels < 10*pi/4*diam^2). Automatically sets MinNpix and MaxNpix below (1/4). 

MaxNpix --- maximum number of pixels per ROI. Automatically set by diameter, if provided. 

MinNpix --- minimum number of pixels per ROI. Automatically set by diameter, if provided. 

Compact --- a compactness criterion for how close pixels are to the center of the ROI. 1 is the lowest possible value, achieved by perfect disks. Best to leave this to a high value (i.e. 2) before the manual sorting stage. 

parent --- these are criteria imposed on the parent cluster (before separating connected regions). These are not currently used during the automated step, but are available in the GUI.

parent.minPixRelVar --- significant regions need to have at least >1/10 the mean variance of all regions

parent.MaxRegions --- if there are more non-significant regions than this number, this parent ROI is probably very spread out over many small components and its connected regions are not good cells: it will be discarded. 

### IX. Example database entry ###

The following is a typical database entry in the local make_db file, which can be modelled after make_db_example. The folder structure assumed is RootStorage/mouse_name/date/expts(k) for all entries in expts(k). 

i = i+1;
db(i).mouse_name    = 'M150329_MP009'; 
db(i).date          = '2015-04-27';
db(i).expts         = [5 6]; % which experiments to process together

Other (hidden) options are described in make_db_example.m, and at the top of run_pipeline.m (set to reasonable defaults), and get_signals_and_neuropil.m (neuropil "surround" option). 

