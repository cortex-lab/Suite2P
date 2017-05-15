function [consistent_rois] = find_consistent_rois(roi_cell)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

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

