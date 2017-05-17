function [overlaps] = cat_overlap(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%% concatenate overlap %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allows daisy-chaining of multiple registers2p regi variables. NB all must
% have been registered to the same time-point 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% -------
% - none: opens a file selection gui, select multiple regi variables

% - cell-array of strings: each cell contains the full file-path to a regi
% variable

% - structure: multiple regi variables saved into one structure with
% fields .rois, .targets and .dat (and their corresponding sub-fields)

% Output: 
% -------
% - overlaps.fnames: cell array where each cell contains the name of an
% F...proc.mat file containing the ROIs from which the correspondingly
% indexed columns in the below fields are drawn (i.e. overlaps.fnames{1}
% gives the filename for overlaps.rois.idcs(:,1))

% - overlaps.rois.idcs: n * m matrix where n = number of overlapping ROIs
% and m = timepoints. Each element of a given column can be directly
% indexed into the variables Fcell & FcellNeu of the F...proc.mat file
% dictated by the corresponding cell in overlaps.fnames (see above).

% - overlaps.rois.iscell_idcs: same as overlaps.rois.idcs except these
% indices correspond to the non-zero elements of the iscell classifier. To
% index into Fcell & FcellNeu do:
%       nn_idcs = find(iscell);
%       Fcell(nn_idcs(overlaps.rois.iscell_idcs(:,1)))

% - overlaps.targets.idcs: n * m matrix where n = number of targets user
% imported to registers2p and m = timepoints. Elements are either NaN if
% target didn't overlap with a suite2P ROI, or else the index of the
% overlapping ROI in Fcell & Fcell Neu of F...proc.mat file dictated by the
% corresponding cell in overlaps.fnames (see above).

% - overlaps.targets.iscell_idcs: Same as overlaps.targets.idcs except
% non-NaN elements are indices of non-zero elements of the iscell
% classifier (see overlaps.rois.iscell_idcs for an example of how to index
% this into Fcell & FcellNeu).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
overlaps = [];
f = {};
d = {};
if ~isempty(varargin)
    if iscell(varargin{1})
        if ischar(varargin{1}{1})
            in = varargin{1};
        end
        for i = 1:numel(in)
            [d{i},ff,ee] = fileparts(in{i});
            f{i} = [ff ee];
        end
    elseif isstruct(varargin{1})
        f = varargin{1};
    end
else
    [f,in] = uigetfile('MultiSelect','on');
    if ~iscell(f)
        f = {f};
    end
    for i = 1:numel(f)
        d{i} = in;
    end
end
if iscell(f)
    rois_idcs = cell(numel(f),1);
    rois_iscell_idcs = rois_idcs;
    targets_idcs = rois_idcs;
    targets_iscell_idcs = rois_idcs;
    overlaps.fname = cell(numel(f)+1,1);
    for i = 1:numel(f)
        load([d{i} filesep f{i}])
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
elseif isstruct(f)
    rois_idcs = cell(numel(f),1);
    rois_iscell_idcs = rois_idcs;
    targets_idcs = rois_idcs;
    targets_iscell_idcs = rois_idcs;
    overlaps.fname = cell(numel(f)+1,1);
    for i = 1:numel(f)
        rois_idcs{i} = f(i).rois.idcs;
        rois_iscell_idcs{i} = f(i).rois.iscell_idcs;
        targets_idcs{i} = f(i).targets.idcs;
        targets_iscell_idcs{i} = f(i).targets.iscell_idcs;
        if i == 1
            [~,overlaps.fname{1}] = fileparts(f(i).dat.files(1).filepath);
            [~,overlaps.fname{2}] = fileparts(f(i).dat.files(2).filepath);
        else
            [~,overlaps.fname{i+1}] = fileparts(f(i).dat.files(2).filepath);
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

