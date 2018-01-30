# Suite2p: fast, accurate and complete two-photon pipeline

Registration, cell detection, spike extraction and visualization GUI. 

Algorithmic details in [http://biorxiv.org/content/early/2016/06/30/061507](http://biorxiv.org/content/early/2016/06/30/061507).

[![IMG](https://img.youtube.com/vi/xr-flH2Ow2Y/0.jpg)](https://www.youtube.com/watch?v=xr-flH2Ow2Y)

This code was written by Marius Pachitariu and members of the cortexlab (Kenneth Harris and Matteo Carandini). It is provided here with no warranty. For support, please open an issue directly on github. 

# Examples
An example dataset (with master_file and make_db) is provided [here](https://drive.google.com/open?id=0B649boZqpYG1amlyX015SG12VU0).

# I. Introduction

This is a complete, automated pipeline for processing two-photon Calcium imaging recordings. It is simple, fast and yields a large set of active ROIs. A GUI further provides point-and-click capabilities for refining the results in minutes. The pipeline includes the following steps

1. **X-Y subpixel registration**: using a modification of the phase correlation algorithm and subpixel translation in the FFT domain. If a GPU is available, this completes in 20 minutes per 1h of recordings at 30Hz and 512x512 resolution.

1. **SVD decomposition**: this provides the input to cell detection and accelerates the algorithm. 

1. **Cell detection**: using clustering methods in a low-dimensional space. The clustering algorithm provides a positive mask for each ROI identified, and allows for overlaps between masks. There is also an option to perform automated red cell detection.

1. **Signal extraction**: by default, all overlapping pixels are discarded when computing the signal inside each ROI, to avoid using "demixing" approaches, which can be biased. The neuropil signal is also computed independently for each ROI, as a weighted pixel average, pooling from a large area around each ROI, but excluding all pixels assigned to ROIs during cell detection. 

1. **Spike deconvolution**: cell and neuropil traces are further processed to obtain an estimate of spike times and spike "amplitudes". The amplitudes are proportional to the number of spikes in a burst/bin. Even under low SNR conditions, where transients might be hard to identify, the deconvolution is still useful for temporally-localizing responses. The cell traces are baselined using the minimum of the (overly) smoothed trace. 

1. **Neuropil subtraction**: coefficient is estimated iteratively together with spike deconvolution to minimize the residual of spike deconvolution. The user is encouraged to also try varying this coefficient, to make sure that any scientific results do not depend crucially on it. 

1. **Automatic and manual curation**: the output of the cell detection algorithm can be visualized and further refined using the included GUI. The GUI is designed to make cell sorting a fun and enjoyable experience. It also includes an automatic classifier that gradually refines itself based on the manual labelling provided by the user. This allows the automated classifier to adapt for different types of data, acquired under different conditions. (*README FOR GUI AT https://github.com/cortex-lab/Suite2P/blob/master/gui2P/README.md*)


# II. Getting started

The toolbox runs in Matlab and currently only supports tiff file inputs. To begin using the toolbox, you will need to make local copies (in a separate folder) of two included files: master_file and make_db. It is important that you make local copies of these files, otherwise updating the repository will overwrite them (and you can lose your files). The make_db file assembles a database of experiments that you would like to be processed in batch. It also adds per-session specific information that the algorithm requires such as the number of imaged planes and channels. The master_file sets general processing options that are applied to all sessions included in make_db, UNLESS the option is over-ridden in the make_db file. 

### Example database entry
Look into make_db_example for more detailed examples.

The following is a typical database entry in the local make_db file, which can be modelled after make_db_example. The folder structure assumed is RootStorage/mouse_name/date/expts(k) for all entries in expts(k).
```
i = i+1; 
db(i).mouse_name = 'M150329_MP009'; 
db(i).date = '2015-04-27'; 
db(i).expts = [5 6]; % which experiments to process together
```
Other (hidden) options are described in make_db_example.m, and at the top of run_pipeline.m (set to reasonable defaults).

### Running the pipeline
Change paths in master_file to the paths to your local toolbox and to your data. Then run this function. The master_file creates the ops0 variable and the db0 variable, and runs the main pipeline:

```
run_pipeline(db, ops);
```

### Spike deconvolution

For spike deconvolution, you need to download the OASIS github (https://github.com/zhoupc/OASIS_matlab) and add the path to this folder on your computer to the top of your master_file 
```
addpath(genpath('pathtoOASIS')))
```
To run spike deconvolution (after running the pipeline), run
```
add_deconvolution(ops0, db);
```

For L0 spike deconvolution, you need to run mex -largeArrayDims SpikeDetection/deconvL0.c (or .cpp under Linux/Mac). If you're on Windows, you will need to install Visual Studio Community in order to mex files in matlab. To choose this deconvolution method, set  
```
ops0.deconvType = 'L0';
```
See this paper comparing spike deconvolution methods for more information on choosing deconvolution methods/parameters: http://www.biorxiv.org/content/early/2017/06/27/156786

You can also run spike deconvolution without running the entire pipeline by calling wrapperDECONV(ops,F,N), where F and N are the fluorescence and neuropil traces respectively, while ops specifies some deconvolution parameters like sampling rate and sensory decay timescale. See the function help for more information.

----------
Below we describe the outputs of the pipeline first, and then describe the options for setting it up, and customizing it. Importantly, almost all options have pre-specified defaults. Any options specified in `master_file` in ops0 overrides these defaults. Furthermore, any option specified in the `make_db` file (experiment specific) overrides both the defaults and the options set in master_file. This allows for flexibility in processing different experiments with different options. The only critical option that you'll need to set is ops0.diameter, or db(N).diameter. This gives the algorithm the scale of the recording, and the size of ROIs you are trying to extract. We recommend as a first run to try the pipeline after setting the diameter option. Depending on the results, you can come back and try changing some of the other options.  

Note: some of the options are not specified in either the example `master_file` or the example make_db file. These are usually more specialized features.

# III. Outputs.

The output is a struct called dat which is saved into a mat file in ResultsSavePath using the same subfolder structure, under a name formatted like `F_M150329_MP009_2015-04-29_plane1`. It contains all the information collected throughout the processing, and contains the ROI and neuropil traces in Fcell and FcellNeu, and whether each ROI j is a cell or not in stat(j).iscell. stat(j) contains information about each ROI j and can be used to recover the corresponding pixels for each ROI in stat.ipix. The centroid of the ROI is specified in stat as well. Here is a summary of where the important results are:

cell traces are in `Fcell`  
neuropil traces are in `FcellNeu`

deconvolved traces are in `sp`

Each cell of the above structures is a different experiment from db.expts.
manual, GUI overwritten "iscell" labels are in `stat.iscell`  
 
stat(icell) contains all other information:  
* _iscell_: automated label, based on anatomy  
* _neuropilCoefficient_: multiplicative coefficient on the neuropil signal
* _xpix, ypix_: x and y indices of pixels belonging to this max. These index into the valid part of the image (defined by ops.yrange, ops.xrange).   
* _ipix_: linearized indices ((ypix, xpix) --> ypix + (xpix-1) x Ly) of pixels belonging to this mask.   
* _isoverlap_: whether the pixels overlap with other masks.     
* _lam, lambda_: mask coefficients for the corresponding pixels. lambda is the same as lam, but normalized to 1.   
* _med_: median of y and x pixels in the ROI (indices for the valid part of the image, defined by ops.yrange, ops.xrange).   
* _blockstarts_: the cumulative number of frames per block. Clould be useful for concatenating experiments correctly (some planes will have fewer frames/block). 
* _footprint, mrs, mrs0, cmpct, aspec_ratio, ellipse, mimgProj, skew, std, maxMinusMed, top5pcMinusMed_: these are used by the automated classifier to label an ROI as cell or not. see section IX for details.

There are fields for red cell detection too (see the section on **Identifying red cells**)

The settings for the registration and the mean image are also output in the `ops` structure:
* _mimg_: target mean image computed at the beginning of registration to which all frames are aligned
* _mimg1_: mean image computed from all the frames across all experiments
* _DS_: offsets computed in XY

# IV. Input-output file paths

* _RootStorage_ --- the root location where the raw tiff files are  stored.
* _RegFileRoot_ --- location on local disk where to keep the registered movies in binary format. This will be loaded several times so it should ideally be an SSD drive. (to view registered movies use script "view_registered_binaries.m" in main folder)
* _ResultsSavePath_ --- where to save the final results. 
* _DeleteBin_ --- deletes the binary file created to store the registered movies
* _RegFileTiffLocation_ --- where to save registered tiffs (if empty, does not save)
**if you want to save red tiffs, then specify ops.RegFileTiffLocation and set ops.REDbinary = 1**

All of these filepaths are completed with separate subfolders per animal and experiment, specified in the make_db file. Your data should be stored under a file tree of the form

\RootStorage\mouse_name\session\block\*.tif(f)

If you don't want to use this folder structure, see the make_db_example file for alternatives. The make_db_example file also shows how to group together tiffs from different experiments (i.e. different subfolders within this folder structure).

The output is a struct called dat which is saved into a mat file in ResultsSavePath using the same subfolder structure, under a name formatted like F_M150329_MP009_2015-04-29_plane1. It contains all the information collected throughout the processing, and contains the fluorescence traces in dat.Fcell and whether a given ROI is a cell or not in dat.stat(N).iscell. dat.stat contains information about each ROI and can be used to recover the corresponding pixels for each ROI N in dat.stat(N).ipix. The centroid of the ROI N is specified in dat.stat(N) as well.

# V. Options

### Registration

* _showTargetRegistration_ --- whether to show an image of the target frame immediately after it is computed. 
* _PhaseCorrelation_ --- whether to use phase correlation (default is phase-correlation, if 0 then cross-correlation used).
* _SubPixel_ --- accuracy level of subpixel registration (default is 10 = 0.1 pixel accuracy)
* _kriging_ --- compute shifts using kernel regression with a gaussian kernel of width 1 onto a grid of 1/SubPixel (default is kriging = 1)
* _maxregshift_ --- maximum amount of movement allowed in FOV (default is 10% of max(y pixels, x pixels))
* _maskSlope_ --- slope on the taper mask applied to image (default is 2 pixel exponential decay)
* _NimgFirstRegistration_ --- number of randomly sampled images to compute target image from
* _NiterPrealign_ --- number of iterations for the target computation (iterative re-alignment of subset of frames)
* _smooth_time_space_ --- convolves raw movie with a Gaussian of specified size in specified dimensions;
                      [t]: convolve in time with gauss. of std t, [t s]: convolve in time and space,
                      [t x y]: convolve in time, and in space with an ellipse rather than circle
* _nimgbegend_ --- number of frames over which to average at beginning and end of experiment (if worried about drift) (default is 0 frames)
                      
**Block Registration (for high zoom/npixels)**

* _nonrigid_ --- set to 1 for non-rigid registration (or set numBlocks > 1)
* _numBlocks_ --- 1x2 array denoting the number of blocks to divide image in y and x (default is [8 1]) 
* _blockFrac_ --- percent of image to use per block (default is 1/(numBlocks-1))          
* _quadBlocks_ --- interpolate block shifts to single line shifts (6 blocks -> 512 lines) by fitting a quadratic function (default is 1)
* _smoothBlocks_ --- if quadBlocks = 0, then smoothBlocks is the standard deviation of the gaussian smoothing kernel

**Bidirectional scanning issues (frilly cells - default is to correct)**

* _dobidi_ --- compute bidirectional phase offset from images (default 1)
* _BiDiPhase_ --- value of bidirectional phase offset to use for computation (will not compute bidirectional phase offset)

**Recordings with red channel**

* _AlignToRedChannel_ --- perform registration to red channel (non-functional channel) rather than green channel
**(this assumes that you have a red channel for all recordings in db.expts)**
* _redMeanImg_ --- compute mean image of red channel from experiments with red and green channel
(you do not need to have a red channel for all db.expts, computes only from db.expred -- so make sure this isn't empty!!)
* _REDbinary_ --- compute a binary file of the red channel (like the green channel binary) from db.expred

output ops.mimgRED will contain mean image (if AlignToRedChannel, redMeanImg or REDbinary = 1)

**Splitting large tiffs for registration if running out of memory (e.g. 2048 x 2048 pixel images)**
currently only works with rigid registration, where each section of FOV is registered separately
* splitFOV --- 1x2 array specifying chunk size in y and x (default is [1 1])

### Cell detection

* _sig_ --- spatial smoothing constant: smooths out the SVDs spatially. Indirectly forces ROIs to be more round. 
* _nSVDforROI_ --- how many SVD components to keep for clustering. Usually ~ the number of expected cells. 
* _ShowCellMap_ --- whether to show the clustering results as an image every 10 iterations of the clustering
* _getROIs_ --- whether to run the ROI detection algorithm after registration
* _stopSourcery_ --- stop clustering if # of ROIs extracted < (ROIs extracted on iteration 1) x (stopSourcery) (default is 1/10)
* _maxIterRoiDetection_ --- maximum number of clustering iterations (default is 100)

### SVD decomposition

* _NavgFramesSVD_ --- for SVD, data has to be temporally binned. This number specifies the final number of points to be obtained after binning. In other words, datasets with many timepoints are binned in higher windows while small datasets are binned less. 
* _getSVDcomps_ --- whether to obtain and save to disk SVD components of the registered movies. Useful for pixel-level analysis and for checking the quality of the registration (residual motion will show up as SVD components). This is a separate SVD decomposition from that done for cell clustering (does not remove a running baseline of each pixel). 
* _nSVDforROI_ --- how many SVD components to keep.
* _getSVDcomps_ --- whether to obtain and save to disk SVD components of the registered movies. Useful for pixel-level analysis and for checking the quality of the registration (residual motion will show up as SVD components). This is a separate SVD decomposition from that done for cell clustering (does not remove a running baseline of each pixel). 
* _nSVD_ --- if _getSVDcomps_=1, then _nSVD_ specifies how many components to save


### Signal extraction

* _signalExtraction_ --- how should the fluorescence be extracted? The 'raw' option restricts cells to be non-overlapping, 'regression' option allows cell overlaps. The neuropil model is a set of spatial basis functions that tile the FOV. The 'surround' option means that the cell's activity is the weighted sum of the detected pixels (weighted by lambda). The neuropil is computed as the sum of activity of surrounding pixels (excluding other cells in the computation).

### Neuropil options

* _ratioNeuropil_ --- used for both spatial basis functions and surround neuropil - the spatial extent of the neuropil as a factor times the radius of the cells (ops.ratioNeuropil * cell radius = neuropil radius)

if using surround neuropil (_signalExtraction_ = 'surround')
* _innerNeuropil_ --- padding in pixels around cell to exclude from neuropil
* _outerNeuropil_ --- radius of neuropil surround (set to Inf to use ops.ratioNeuropil)
* _minNeuropilPixels_ --- minimum number of pixels in neuropil surround (radius will expand until this number is reached)

### Spike deconvolution 

* _imageRate_ --- imaging rate per plane. 
* _sensorTau_ --- decay timescale (in seconds).
* _maxNeurop_ --- neuropil contamination coef has to be less than this (sometimes good to impose a ceiling at 1, i.e. for interneurons)
* _deconvType_ --- which type of deconvolution to use (either 'L0' or 'OASIS') 

### Identifying red cells

use function `identify_redcells_sourcery(db, ops0)` to identify cells with red

Outputs are appended to stat in F.mat file
* redratio = red pixels inside / red pixels outside
* redcell = redratio > mean(redratio) + _redthres_ x std(redratio)
* notred = redratio < mean(redratio) + _redmax_ x std(redratio)

Options are 
* _redthres_  --- higher thres means fewer red cells (default 1.35)
* _redmax_ --- the higher the max the more NON-red cells (default 1)

### Measures used by classifier

The Suite2p classifier uses a number of features of each ROI to assign cell labels to ROIs. The classifier uses a naive Bayes approach for each feature, and models the distribution of each feature with a non-parametric, adaptively binned empirical distribution. The classifier is initialized with some standard distributions for these features, but is updated continuously with new data samples as the user refines the output manually in the GUI. 

The features used are the following (can see values for each ROI by selecting it in the GUI). 
* std = standard deviation of the cell trace, normalized to the size of the neuropil trace  
* skew = skewness of the neuropil-subtracted cell trace  
* pct = mean distance of pixels from ROI center, normalized to the same measuree for a perfect disk  
* footprint = spatial extent of correlation between ROI trace and nearby pixels  
* mimgProjAbs = whether this ROI shape is correlated to the shape on the mean image  
* aspect_ratio = of an ellipse fit to the ROI  


