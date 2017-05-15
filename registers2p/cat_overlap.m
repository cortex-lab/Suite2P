function [overlaps] = cat_overlap(varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

overlaps = [];
f = [];
if ~isempty(varargin)
    if iscell(varargin{1})
        if ischar(varargin{1}{1})
            f = varargin{1};
        end
    else
        if ischar(varargin{1}{1})
            f{1} = varargin{1};
        end
    end
else
    [f,d] = uigetfile('MultiSelect','on');
    if ~iscell(f)
        f = {f};
    end
end
if ~isempty(f)
    rois_idcs = cell(numel(f),1);
    rois_iscell_idcs = rois_idcs;
    targets_idcs = rois_idcs;
    targets_iscell_idcs = rois_idcs;
    overlaps.fname = cell(numel(f)+1,1);
    for i = 1:numel(f)
        load([d filesep f{i}])
        rois_idcs{i} = regi.rois.idcs;
        rois_iscell_idcs{i} = regi.rois.iscell_idcs;
        targets_idcs{i} = regi.targets.idcs;
        targets_iscell_idcs{i} = regi.targets.iscell_idcs;
        if i == 1
            [~,overlaps.fname{1}] = fileparts(regi.dat.files(1).filepath);
            [~,overlaps.fname{2}] = fileparts(regi.dat.files(2).filepath);
        else
            [~,overlaps.fname{i+1}] = fileparts(regi.dat.files(2).filepath);
        end
    end
    overlaps.rois.idcs = find_consistent_rois(rois_idcs);
    overlaps.rois.iscell_idcs = find_consistent_rois(rois_iscell_idcs);
    if max(cellfun('isempty',targets_idcs)) == 0
        overlaps.targets.idcs = [targets_idcs{1}(:,1) cell2mat(cellfun(@(x) x(:,2),targets_idcs,'UniformOutput',0)')];
        overlaps.targets.iscell_idcs = [targets_iscell_idcs{1}(:,1) cell2mat(cellfun(@(x) x(:,2),targets_iscell_idcs,'UniformOutput',0)')];
    else
        overlaps.targets.idcs = [];
        overlaps.targets.iscell_idcs = [];
    end
else
    warning('Import error - either provide no input and use file-selection gui or provide a cell-array of filepath strings')
end
end

