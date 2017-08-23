# How to use

1) Run new_main from the command line in the gui2P folder. 

2) Load single to load a single result file. Find a file formatted as F_M150329_MP009_2015-04-29_plane1_Nk650 in your results path. 

3) After <20 seconds this will load the results, providing you with a left image of ROIs classified as cells, and a right window of non-classified ROIs. The color is randomly chosen for each ROI.

The most important capability of the GUI is the ability to point and click on a pixel in either of the two views to select its cluster. The fluorescence trace will then be displayed. Right-clicking on an ROI in the left window pushes that ROI to the right window and viceversa. This allows quick manual sorting of the detected clusters and is the main function provided by this GUI. For RED cell classification, use the middle click (the mouse wheel) to change a cell from RED to NOTRED and viceversa.

4) Once you are done clicking on cells, press "Save"  to save the results back into a new file, identical to the file you have loaded from with a '_proc' at the end of the file name, and with a new field iscell in the stat structure, which indicates the result of the GUI. A few other settings of the GUI are saved, so that the session can be resumed at a later time. 

## there are help buttons in the GUI too (DON'T PANIC and HELP buttons in red) 

# GUI FEATURES
## Classifier

The classifier builds a non-parametric model of the cell/non-cell parameters (compactness, skewness, size etc.). The model is updated after every user manual session (after you press save). You can have different classifiers for different types of data (i.e. somas, dendrites, boutons, different brain areas, calcium indicators, or zoom levels). The classifier assigns the initial "iscell" labels to the ROIs. Initially, it might not work very well, but will improve as you make choices in the GUI (and save the "proc" files). 

How to use it: 
To begin, press "new classifier?" and choose one of the provided priors.
The last used classifier will be automatically selected, every time you load a new plane.
As you process datasets, a database of your cells is built, and the classifier is re-trained.
The prior counts for about 300 cells, so your choices will start making a difference after a few hundred manually validated cells.


## Background

ROI: shows the masks identified by Suite2p, colored according to the property selected under "mask color"

CORR: shows the correlation map between a pixel and its nearby pixels

MEAN: average registered image

RED: average registered image of "red"/secondary color channel

RED - GREEN: subtracts off the contamination from "green"/primary channel

PROC: switch to toggle image contrast normalization

   hint: the letters in paranthesis are keyboard shortcuts


## Mask Colors

This applies only if you select "ROIs" under "background"

RANDOM: color chosen randomly - RED cells labelled in red (if secondary channel)

   *** middle mouse-click on an ROI to switch RED label ON/OFF ***

for all other selections, color/hue for each ROI varies purple to yellow (scale bar for colors shown below buttons)

CLASSIFIER: probability assigned by classifier

SKEW: skewness of activity, after neuropil correction and some smoothing

MEANIMG: weighting of activity mask onto mean image

CMPCT: compactness of ROI pixels. Smallest is 1, for disks

FOOT: "footprint" of ROI; ~ number of correlated neighboring pixels

RED: probability of being a red-tagged cell, assigned by algorithm in identify_redcells_sourcery.m

   hint: the letters in paranthesis are keyboard shortcuts

## Mask Brightness

The Mask Brightness can be used to switch between different ways of displaying the masks. The algorithm returns continuous valued masks that are normalized to 1 per ROI. 

MaskBrightness= UnitVectors encodes the mask weights into the brightness/value of the image and is the default option. 

MaskBrightness = VarianceExp uses the global variance explained of each pixel by its ROIs. This setting can be used to see the global magnitude of variance in each ROI. Bright ROIs will have a lot of variance and transients. Importantly, the variance is computed on high-pass filtered data, so ROIs with drifting fluorescence will not appear to have a lot of variance. 

## ZOOM

The quadrant view allows you to zoom into pre-specified 3 by 3 quadrants of the image. All the options of the GUI are available in the zoomed view. 


