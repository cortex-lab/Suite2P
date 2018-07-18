function [h]= identify_classifier(h)
% load classifier files

rootS2p = which('run_pipeline');
if isempty(rootS2p)
    warndlg('Suite2p location not in PATH! where is run_pipeline.m?')
    [filename1,rootS2p]=uigetfile(root, 'Select Data File');
    
end
rootS2p = fileparts(rootS2p);
rootS2p = fullfile(rootS2p, 'configFiles');

flag=0
if isfield(h.dat.ops, 'classifier') && ~isempty(h.dat.ops.classifier)
    if ~exist(fullfile(rootS2p, h.dat.ops.classifier), 'file')
         warndlg('specified classifier database not found, reverting  to last used');
    else
        def_file = h.dat.ops.classifier;
        flag = 1;
        h.dat.cl.fpath  = fullfile(rootS2p, def_file);
        
        hload = load(h.dat.cl.fpath);
        if ~isfield(hload, 'st') || ~isfield(hload, 'statLabels') || ~isfield(hload, 'prior')
            error('found a non-classifier file in configFiles, called %s. \n Please remove and try again!', def_file)
        end

        h.st        = hload.st;
        h.prior     = hload.prior;
        h.statLabels = hload.statLabels;
    end
else
%     warning('no specified classifier database, reverting  to last used');
end
if (flag==0)
    fs = dir(fullfile(rootS2p, '*.mat'));
    if isempty(fs)
        warndlg('no classifier found, please make a new one!')
        h           = make_classifier(h);
    else
        [~, imax]    = max([fs.datenum]);
        
        def_file = fs(imax).name;

        h.dat.cl.fpath  = fullfile(rootS2p, def_file);
        
        hload = load(h.dat.cl.fpath);
        if ~isfield(hload, 'st') || ~isfield(hload, 'statLabels') || ~isfield(hload, 'prior')
            error('found a non-classifier file in configFiles, called %s. \n Please remove and try again!', def_file)
        end

        h.st        = hload.st;
        h.prior     = hload.prior;
        h.statLabels = hload.statLabels;
        h.dat.ops.classifier = def_file;
    end
end

end

function h = make_classifier(h)
% new classifier button
rootS2p = which('run_pipeline');
if isempty(rootS2p)
    warndlg('Suite2p location not in PATH! where is run_pipeline.m?')
    [filename1,rootS2p]=uigetfile(root, 'Select run_pipeline.m');
    
end
rootS2p = fileparts(rootS2p);
rootS2p = fullfile(rootS2p, 'configFiles');

def_name = fullfile(rootS2p, 'cl_new.mat');
[FileName,PathName] = uiputfile('*.mat', 'Create new classifier', def_name); 

if FileName
    load(fullfile(rootS2p, 'priors', 'prior_default.mat'));
    st = [];
    save(fullfile(PathName, FileName), 'st', 'prior', 'statLabels')
    
    h.st = st;
    h.prior = prior;
    h.statLabels = statLabels;
    
    h.dat.cl.fpath          = fullfile(PathName, FileName);
    h                       = classROI(h);
    
%     h = splitROIleftright(h);
%     h = buildLambdaValue(h);
%     redraw_figure(h);
    
    hload = load(h.dat.cl.fpath);
    h.st        = hload.st;
    h.statLabels = hload.statLabels;

    h.prior     = hload.prior;
else
    error('you must make a new classifier first!')
end

end
