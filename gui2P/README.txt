1) Run new_main from the command line in the gui2P folder. 

2) Load single to load a single result file. Find a file formatted as F_M150329_MP009_2015-04-29_plane1_Nk650 in your results path. 

3) After <20 seconds this will load the results, providing you with a left image of ROIs classified as cells, and a right window of non-classified ROIs. The color is randomly chosen for each ROI.

4) The most important capability of the GUI is the ability to point and click on a pixel in either of the two views to select its cluster. The fluorescence trace will then be displayed. Right-clicking on an ROI in the left window pushes that ROI to the right window and viceversa. This allows quick manual sorting of the detected clusters and is the main function provided by this GUI. 

5) The Mask Type can be used to switch between different ways of displaying the masks. The algorithm returns continuous valued masks that are normalized to 1 per ROI. MaskType= UnitVectors encodes the mask weights into the brightness/value of the image and is the default option. MaskType = VarianceExp uses the global variance explained of each pixel by its ROIs. This setting can be used to see the global magnitude of variance in each ROI. Bright ROIs will have a lot of variance and transients. Importantly, the variance is computed on high-pass filtered data, so ROIs with drifting fluorescence will not appear to have a lot of variance. 

6) The quadrant view allows you to zoom into pre-specified 3 by 3 quadrants of the image. All the options of the GUI are available in the zoomed view. 

7) The Display view allows you to cycle between the clustering view (ROIs), the correlation map, the mean image of the recording and the mean RED channel, if available. These same views can be selected from the keyboard under buttons Q, W, E,R respectively. 

8) The classifier builds a non-parametric model of the cell/non-cell parameters (compactness, skewness, size etc.). The model is updated after every user manual session (after you press save). You can have different classifiers for different types of data. 

9) Once you are done, press "Save"  to save the results back into a new file, identical to the file you have loaded from with a '_proc' at the end of the file name, and with a new field iscell in the stat structure, which indicates the result of the GUI. A few other settings of the GUI are saved, so that the session can be resumed at a later time. 
