Under construction, check back in 24 hours.....

This is a complete automated pipeline for processing two-photon Calcium imaging recordings. It is very simple, very fast and yields a large set of active ROIs. A GUI further provides point-and-click capabilities for refining the results and obtaining a final set of ROIs in minutes. The pipeline includes the following steps

1) X-Y subpixel registration --- using a version of the phase correlation algorithm and subpixel translation in the FFT domain. If a GPU is available, this completes in 20 minutes per 1h of recordings at 30Hz and 512x512 resolution.

2) SVD decomposition --- this provides the basis for a number of pixel-level visualization tools. 

3) Cell detection --- using clustering methods in a low-dimensional space of the fluorescence activity.

4) Manual curation --- the output of the cell detection algorithm can be visualized and further refined using a GUI available at https://github.com/marius10p/gui2P/. The GUI is designed to make cell sorting a fun and enjoyable experience. 

The toolbox currently only supports tiff file inputs. To begin using the toolbox, you will need to make local copies (in a separate folder) of two 
included files: master_file and make_db. The make_db file assembles a database of experiments that you would like to be processed in batch. It also 
adds per-session specific information that the algorithm requires such as the number of imaged planes and channels. The master_file sets general 
processing options that will be applied to all sessions included in make_db, UNLESS the option is over-ridden in the make_db file.  The global and session-specific options are described in detail below. 