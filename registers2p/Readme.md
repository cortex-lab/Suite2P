# Register S2P: semi-automatic cross-day ROI registration #

[![IMG](https://img.youtube.com/vi/AwVYeyh2O6A/0.jpg)](https://youtu.be/AwVYeyh2O6A)

Allows users to import two Suite2P output structures (F...proc.mat files) generated from a single FOV across 2 time-points (e.g. > 1 day) and semi-automatically register the two images of the FOV and the corresponding Suite2P ROIs. It also has the optional functionality of loading in user-defined target centroids which can then be mapped onto Suite2P ROIs (and overlapped across days). Allows quick and easy manual curation of overlaps. Use cat_overlap.m function provided to daisy-chain registered ROIs recorded at multiple time-points which have all been registered to one reference time-point.

---

### 1. Image window (left) ###
- Displays 1 of the 2 imported datasets: Q = dataset 1, W = dataset 2.
- Can display 1 of 5 masks for each dataset: 1 - 5 correspond to mask-types listed about the image window.
- Hit P to toggle processing on any of these images (enhances contrast).
- Can display ROIs of the displayed dataset: A = dataset 1 ROIs, S = dataset 2 ROIs.
- Clicking: Click once to open a lasso tool to select user-imported target centroids (if present). The cursor will switch to a cross-hair. Subsequently click and drag to define the lasso area around targets you want to select. Double-click to finish defining area. Selected targets within the area will turn from red to yellow. Subsequent clicks and drags will move targets around. Hit Enter to de-select targets and leave them where you have dragged them. If you select targets you can revert them to their default locations by pressing D.
---

### 2. ROI window (right) ###
- Displays ROIs of both images overlaid with grey-level dictating which image is currently displayed on the Image window: white = displayed image, grey = non-displayed image. This is also dictated by Q (dataset 1) and W (dataset 2).
- ROIs can be in two states: non-overlapping (open contour) or overlapping (filled contour). All ROIs start in default non-overlapping state until decided otherwise by the user.
- Clicking: Clicking ROIs will toggle the overlap state of selected ROIs that have >0 pixel overlap with each other. In the condition where ROI 1 overlaps with 2 ROIs clicking it will cause it to overlap with the ROI with which it shares the most pixels. An ROI can only overlap with max 1 other ROI so designating ROI 1 (which already overlaps with ROI 2) to overlap with ROI 3 will result in ROI 2 to revert to the non-overlapping state.
---

### 3. Registration ###
- Import 2 Suite2P F...proc.mat files  (Load dataset 1, Load dataset 2).
- Click Register.
- Navigate the window that opens, clicking corresponding control points in the left and right image (in order, i.e. click point 1 in left-hand image, click 1 in right-hand image, click 2 in left-hand image etc.).
- Click "File -> Close Control Point Selection Tool".
- Wait for transform to be applied.
- Note that whichever mask and processing state you are currently viewing will be used when registering. This cannot be changed during the registering process.
---

### 4. Detect ROI overlap ###
- In "ROIs" panel above image window 2 decide on proportion of pixel overlap you require (this is calculated from the ROIs in dataset 1, i.e. if an ROI from dataset 2 covers all pixels in the ROI from dataset 1 then this corresponds to 100% overlap).
- Click "Detect overlap".
- Manually curate overlap states of ROIs by clicking them (see sections 1 & 2) in ROI display (right). NB can use O to toggle cursor mirroring on the two displays (actual cursor position on either display will also be displayed on the other display as a blue circle).
- Clicking "Reset" will reset all ROI overlap states to non-overlapping.
---

### 5. Detect targets overlap ###
- Save a target centroids file in one of 2 formats: a black (zeros) .tiff/.tif file of the same dimensions as the Suite2P datasets with non-zero elements at the pixel locations of your targets, or a .mat file containing a variable called "params" with a field "targets_yxzc" which is an n * 2 matrix where n is the number of targets, column 1 is the y co-ordinates and column 2 is the x co-ordinates of each target.
- Targets will be plotted as small red dots surrounded by one small and one large concentric circle. Concentric circles are stationary and indicate default positions of targets as imported.
- In the "Targets" panel above image window 2 decide on the minimum pixel distance you require between each of your targets and surrounding ROI centroids.
- Click "Find closest ROIs".
- Manually curate overlap of targets with ROIs by clicking, lassoing around desired targets and dragging/nudging (see section 1).
- Selected targets will become yellow. Targets that overlap with an ROI on one or other day will cause the border of that ROI to toggle to red. Targets the overlap with ROIs on both days will cause both ROIs' borders to toggle to red and the target marker will increase in size.
- Clicking "Reset" will reset all targets to be non-overlapping with ROIs.
---

### 6. Outputs ###
Click "Save analysis..." to save registered ROIs (and targets) for subsequent analysis in a variable called regi with fields:
- rois.idcs = n * 2 matrix where n is the number of overlapping rois and columns correspond to the 2 time-points. Elements are the indices of raw Suite2P ROIs in the F...proc.mat file, not parsed by the iscell classifer. I.e. an element that equals 10 in rois.idcs is truly ROI 10 in F...proc.mat. Use this to index directly into dat.Fcell/dat.FcellNeu.
- rois.iscell_idcs = n * 2 matrix where n is the number of overlapping rois and columns correspond to the 2 time-points. Elements are the indices of Suite2P ROIs in the F...proc.mat file after being parsed by the iscell classifier. If you are only importing positively classified ROIs (iscell == 1) then use this to index into your resulting matrix directly: classfied_cells(rois.iscell_idcs,:). If you are importing all ROIs (iscell==1 | iscell == 0) then index into this matrix using the non-negative indices of iscell: nn_idcs = find(iscell); nonclassified_cells(nn_idcs(rois.iscell_idcs),:).
- targets.idcs = n * 3 matrix where n is the number of targets that you imported, columns 1 and 2 are the two time-points and column 3 is 0 | 1 depending on whether the target has overlapping ROIs across the 2 time-points (1) or only has an ROI on one/other/neither time-point (0). Elements in the first 2 columns are NaN if no ROI overlapped with that target, or indices of raw Suite2P ROIs not parsed by the iscell classifier (see roi.idcs above) if an ROI did.
- targets.iscell_idcs = n * 3 matrix where n is the number of targets that you imported, columns 1 and 2 are the two time-points and column 3 is 0 | 1 depending on whether the target has overlapping ROIs across the 2 time-points (1) or only has an ROI on one/other/neither time-point (0). Elements in the first 2 columns are NaN if no ROI overlapped with that target, or indices of Suite2P ROIs after they have been parsed with the iscell classifier (see roi.iscell_idcs above) if an ROI did.
---

### 7. Controls ###
- O: Toggle cursor mirror on left and right displays
- Q: Display mean image (left-display) and highlight ROIs (right-display) from dataset 1
- W: Display mean image (left-display) and highlight ROIs (right-display) from dataset 2 
- A: Display ROIs (left-display) and highlight ROIs (right-display) from dataset 1
- S: Display ROIs (left-display) and highlight ROIs (right-display) from dataset 2 
- P: Toggle processing on currently displayed mask (left-display)
- 1 - 5: Displays number mask (listed above image display) of currently displayed dataset (left-display)
- Enter: When user-imported centroid targets are selected this deselects them (leaving them where they currently are)
- D: When user-imported centroid targets are selected this returns them to their default location
- Arrow keys: When user-imported centroid targets are selected this will nudge them around the image
- Left-click (left-display): click once to open a lasso tool to select user-imported target centroids (if present). The cursor will switch to a cross-hair. Subsequently click and drag to define the lasso area around targets you want to select. Double-click to finish defining area. Selected targets within the area will turn from red to yellow. Subsequent clicks and drags will move targets around. Hit Enter to de-select targets and leave them where you have dragged them.
- Left-click (right-display): toggle the overlap state of selected ROIs that have >0 pixel overlap with each other.
---

### 7. Additional functions ###
- cat_overlap: Allows daisy-chaining of multiple registers2p regi variables. NB all must have been registered to the same time-point 1. See function description for details of inputs and outputs.
- find_consistent_rois: sub-function of cat_overlap, essentially similar functionality but takes a cell array constructed from sub-fields of regi structures instead of the regi structures themselves.
