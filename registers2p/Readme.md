# Register S2P: semi-automatic cross-day Suite2P ROI registration #

[![IMG](https://img.youtube.com/vi/6jutIbOM4Lg/0.jpg)](https://www.youtube.com/watch?v=6jutIbOM4Lg)

- Allows users to import two Suite2P output structures (F..._proc.mat files) generated from a single FOV across 2 time-points (e.g. > 1 day) and semi-automatically register the mean image of the FOV and the corresponding Suite2P ROIs. It also has the optional functionality of loading in user-defined target centroids which can then be mapped onto Suite2P ROIs (and overlapped across days).
---------------------------------------------------------------------------------

### 1. Mean image window (left) ###
- Displays 1 of the 2 imported datasets: Q = dataset 1, W = dataset 2
- Can display ROIs of the displayed dataset: A = dataset 1 ROIs, S = dataset 2 ROIs
- Clicking: Click once to open a lasso tool to select user-imported target centroids (if present). The cursor will switch to a cross-hair. Subsequently click and drag to define the lasso area around targets you want to select. Double-click to finish defining area. Selected targets within the area will turn from red to yellow. Subsequent clicks and drags will move targets around. Hit Enter to de-select targets and leave them where you have dragged them. If you select targets you can revert them to their default locations by pressing D.
---------------------------------------------------------------------------------

### 2. ROI window (right) ###
- Displays ROIs of both images overlaid with Grey-level dictating which image is currently displayed on the Mean image window: white = displayed image, grey = non-displayed image. This is also dictated by Q (dataset 1) and W (dataset 2).
- ROIs can be in two states: non-overlapping (open contour) or overlapping (filled contour). All ROIs start in default non-overlapping state until decided otherwise by the user.
- Clicking: Clicking ROIs will toggle the overlap state of selected ROIs that have >0 pixel overlap with each other. In the condition where ROI 1 overlaps with 2 ROIs clicking it will cause it to overlap with the ROI with which it shares the most pixels. An ROI can only overlap with max 1 other ROI so designating ROI 1 (which already overlaps with ROI 2) to overlap with ROI 3 will result in ROI 2 to revert to the non-overlapping state.
---------------------------------------------------------------------------------

### 3. Registration ###
- Import 2 Suite2P proc.mat files (Load dataset 1, Load dataset 2)
- Click Register
- Navigate the window that opens clicking corresponding control points in the left and right image (in order, i.e. click point 1 in left-hand image, click 1 in right-hand image, click 2 in left-hand image etc.)
- Click "File -> Close Control Point Selection Tool"
- Wait for transform to be applied
---------------------------------------------------------------------------------

### 4. Detect overlap ###
- In "ROIs panel" above image window 2 decide on proportion of pixel overlap you require (this is calculated from the ROIs in dataset 1, i.e. if and ROI from dataset 2 covers all pixels in the ROI from dataset 1 then this corresponds to 100% overlap)
- Click "Detect overlap"
- Manually curate overlap states of ROIs (see sections 1 & 2) in ROI display (right). NB can use O to toggle cursor mirroring on the two displays (cursor position on ROI display will be shown on image display also)
- Clicking "Reset" will reset all ROI overlap states to non-overlapping.
---------------------------------------------------------------------------------

### 5. Controls ###
O: Toggle cursor mirror on left and right displays

Q: Display mean image (left-display) and highlight ROIs (right-display) from dataset 1

W: Display mean image (left-display) and highlight ROIs (right-display) from dataset 2 

A: Display ROIs (left-display) and highlight ROIs (right-display) from dataset 1

S: Display ROIs (left-display) and highlight ROIs (right-display) from dataset 2 

Enter: When user-imported centroid targets are selected this deselects them (leaving them where they currently are)

D: When user-imported centroid targets are selected this returns them to their default location

Arrow keys: When user-imported centroid targets are selected this will nudge them around the image

Left-click (left-display): lick once to open a lasso tool to select user-imported target centroids (if present). The cursor will switch to a cross-hair. Subsequently click and drag to define the lasso area around targets you want to select. Double-click to finish defining area. Selected targets within the area will turn from red to yellow. Subsequent clicks and drags will move targets around. Hit Enter to de-select targets and leave them where you have dragged them.

Left-click (right-display): toggle the overlap state of selected ROIs that have >0 pixel overlap with each other.
