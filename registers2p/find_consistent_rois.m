function [consistent_rois] = find_consistent_rois(roi_cell)

% Used to find consistent overlapping ROIs across > 2 timepoints.
%
% Takes as input a cell-array of length t where cells represents the
% overlap of ROIs from t timepoints compared to a consistent timepoint 1.
% Each cell therefore contains an m * 2 matrix where m = number of
% overlapping ROIs and columns correspond to 2 time-points and each element
% is the index of a given ROI. Time-point 1 (column 1) in all cells must be
% the same.

% This function then finds all time-point 1 ROIs that are consistent across
% all cells (the intersection) and then for each cell finds the ROIs that
% correspond to these initial ROIs. 

% This results in an output matrix "consistent_rois" which is n * (t+1)
% where n is some consistent subset of the union of the first column of
% each cell in the input cell array and t+1 corresponds to the number of
% time-points in the original cell-array + the initial time-point (which
% should be consistent across all cells). Effectively:
%
% consistent_rois(:,1) = overlapping ROIs from timepoint 1
% consistent_rois(:,2) = overlapping ROIs from timepoint 2
%
% consistent_rois(1,:) = ROI 1 across all t time-points

if numel(roi_cell) > 1
    all_initial_rois = cellfun(@(x) {x(:,1)},roi_cell);
    if numel(all_initial_rois) > 1
        consistent_initial_rois = intersect_n(all_initial_rois);
    else
        consistent_initial_rois = all_initial_rois{1};
    end
    consistent_rois = zeros(numel(consistent_initial_rois),numel(all_initial_rois));
    for c = 1:numel(consistent_initial_rois)
        consistent_rois(c,:) = cellfun(@(x) x(x(:,1)==consistent_initial_rois(c),2),roi_cell);
    end
    consistent_rois = [consistent_initial_rois consistent_rois];
else
    consistent_rois = roi_cell{1};
end

end

function [b] = intersect_n(varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if iscell(varargin{1})
    in = varargin{1};
end
a = in{1};
b = in{2};
for i = 3:numel(in)
    b = intersect(a,b,'stable');
    a = in{i};
end
b = intersect(a,b,'stable');

end

