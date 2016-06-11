Under construction, check back in 24 hours.....

Preliminary instructions for a typical cell sorting session.

1) Run new_main from the command line in the gui2P folder. 

2) Load single to load a single result file. Find a file formatted as F_M150329_MP009_2015-04-29_plane1_Nk650 in your results path. Currently the other loading options are unavailable. 

3) After <20 seconds this will load the results, providing you with a left image of ROIs classified as cells, and a right window of non-classified ROIs. The color is randomly chosen for each ROI, and it can be changed by pressing the Randomize Hue button. 

4) The most important capability of the GUI is the ability to point and click on a pixel in either of the two views to select its cluster. The fluorescence trace will then be displayed (ignore the options for processing the fluorescence, they are not implemented yet). Right-clicking on an ROI in the left window pushes that ROI to the right window and viceversa. This allows quick manual sorting of the detected clusters and is the main function provided by this GUI. 

5) The Mask Type can be used to switch between different ways of displaying the masks. The algorithm returns continuous valued masks that are normalized to 1 per ROI. MaskType= UnitVectors encodes the mask weights into the brightness/value of the image and is the default option. MaskType = VarianceExp uses the global variance explained of each pixel by its ROIs. This setting can be used to see the global magnitude of variance in each ROI. Bright ROIs will have a lot of variance and transients. Importantly, the variance is computed on high-pass filtered data, so ROIs with drifting fluorescence will not appear to have a lot of variance. 

6) The quadrant view allows you to zoom into pre-specified 3 by 3 quadrants of the image. All the options of the GUI are available in the zoomed view. 

7) The Display view allows you to cycle between the clustering view (ROIs), the mean image of the recording and the mean RED channel, if available. These same views can be selected from the keyboard under buttons Q, W and E respectively. 

8) ROI rule files are currently unavailable, check back later. 

9) ROI selection tools allow you to globally change the rules for selecting ROIs. Compactness measures how closely the pixels of a ROI are clustered to its centroid. The minimum and best value is 1, achieved only by a perfect disk-shaped ROI. Values less than 1.1-1.2 are typically quite round ROIs. For a 1um pixel, cells have typical pixel counts between 50 and 200.

10) Parent selection rules. Like the ROI selection tools, but the rules are applied to the parent of each ROI. Sometimes, but not always, these rules can further help to classify ROIs as cells. Each ROI is a connected region that comes from a cluster detected during the optimization. This cluster is called the parent cluster, and has typically more than 1 subregions, for example a cell ROI + some noise pixels in a completely different part of the FOV, or two cell ROIs put into the same cluster, or a cell and its dendrites etc. Region count < 40 is a good value and would seem to have the most impact on the classification but you can try the other settings too.

11) Once you are done, press "Save ROIs and F"  to save the results back into the file they have been loaded for. At any future time you can reload this file and work more on the segmentation if you are unhappy with it. 