function varargout = new_main(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @new_main_OpeningFcn, ...
                   'gui_OutputFcn',  @new_main_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before new_main is made visible.
function new_main_OpeningFcn(hObject, eventdata, h, varargin)
% This function has no output args, see OutputFcn.
% varargin   command line arguments to new_main (see VARARGIN)

axes(h.axes2);
set(gca, 'xtick', [], 'ytick', [])
set(gca, 'xcolor', 'w', 'ycolor', 'w')
axes(h.axes3);
set(gca, 'xtick', [], 'ytick', [])
set(gca, 'xcolor', 'w', 'ycolor', 'w')
axes(h.axes4);
set(gca, 'xtick', [], 'ytick', [])
set(gca, 'xcolor', 'w', 'ycolor', 'w')

h.output = hObject;
guidata(hObject, h);
% UIWAIT makes new_main wait for user response (see UIRESUME)
% uiwait(h.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = new_main_OutputFcn(hObject, eventdata, h) 
% varargout  cell array for returning output args (see VARARGOUT);
varargout{1} = h.output;

% main LOAD call
function pushbutton17_Callback(hObject, eventdata, h)
% [filename1,filepath1]=uigetfile('\\zserver\Lab\Share\Marius\', 'Select Data File');

flag = 0;
try
    if isfield(h, 'dat') && isfield(h.dat, 'filename')
        root = fileparts(h.dat.filename);
    else
        root = 'D:\DATA\F\';
    end
    [filename1,filepath1]=uigetfile(root, 'Select Data File');
    h.dat = load(fullfile(filepath1, filename1));
    set(h.figure1, 'Name', filename1);

    flag = 1;
catch
end

if flag
    % if the user selected a file, do all the initializations
rng('default')

% keyboard;
if isfield(h.dat, 'dat')
    h.dat = h.dat.dat;
else
    h.dat.filename = fullfile(filepath1, filename1);
    h.dat.cl.Ly       = numel(h.dat.ops.yrange);
    h.dat.cl.Lx       = numel(h.dat.ops.xrange);
    
    % make up iclut here
    try
        [h.dat.res.iclust, h.dat.res.lambda, h.dat.res.lambda0] =...
            getIclust(h.dat.stat, h.dat.cl);
    catch
    end
    h.dat.res.iclust = reshape(h.dat.res.iclust, h.dat.cl.Ly, h.dat.cl.Lx);
    h.dat.res.lambda = reshape(h.dat.res.lambda, h.dat.cl.Ly, h.dat.cl.Lx);
    
    h.dat.ops.Nk = numel(h.dat.stat);
    h.dat.cl.rands_orig   = .1 + .8 * rand(1, h.dat.ops.Nk);
    h.dat.cl.rands        = h.dat.cl.rands_orig;
    
    if isfield(h.dat.ops, 'clustrules')
       h.dat.clustrules = h.dat.ops.clustrules; 
    end
    
    % set up classifier
    h.dat.cl.threshold  = 0.5;
    h                   = identify_classifier(h);    
    h                   = classROI(h);

    
%     set(h.edit35,'String', num2str(h.dat.res.Mrs_thresh));
%     set(h.edit40,'String', num2str(h.dat.cl.npix_low));
    
    % set all quadrants as not visited
    h.quadvalue = zeros(3);
    for j = 1:3
        for i = 1:3
            set(h.(sprintf('Q%d%d', j,i)), 'BackgroundColor',[.92 .92 .92]);
        end
    end
    
    % make lambda here
    lam = h.dat.res.lambda;
    h.dat.img0.V = max(0, min(1, .75 * ...
        reshape(lam, h.dat.cl.Ly, h.dat.cl.Lx)/mean(lam(lam>1e-10))));
    
    h.dat.ylim = [0 h.dat.cl.Ly];
    h.dat.xlim = [0 h.dat.cl.Lx];    
    
    if ~isfield(h.dat.cl,'redcell')
        h.dat.cl.redcell = zeros(h.dat.ops.Nk, 1);
    end
    h                = splitROIleftright(h);
    
   
    h.dat.F.ichosen = 1; %ceil(rand * numel(icell))
    
    Sat = ones(h.dat.cl.Ly, h.dat.cl.Lx);
    Sat(h.dat.res.iclust==h.dat.F.ichosen) = 0;
    h.dat.img1.Sat     = Sat;
    h.dat.img2.Sat     = Sat;
    
    h = buildHue(h);
    h = buildLambdaValue(h);
    
    % loop through redcells and set h.dat.cl.rands(h.dat.F.ichosen) = 0
    for j = find(h.dat.cl.redcell)
        h.dat.cl.rands(j) = 0;
    end
   
    % x and y limits on subquadrants
    h.dat.figure.x0all = round(linspace(0, 19/20*h.dat.cl.Lx, 4));
    h.dat.figure.y0all = round(linspace(0, 19/20*h.dat.cl.Ly, 4));
    h.dat.figure.x1all = round(linspace(1/20 * h.dat.cl.Lx, h.dat.cl.Lx, 4));
    h.dat.figure.y1all = round(linspace(1/20 * h.dat.cl.Ly, h.dat.cl.Ly, 4));
    
end

% activate all pushbuttons
pb = [84 93 101 86 87 89 90 92 103 98 95 96 102 99 100 1 2];
for j = 1:numel(pb)
    set(eval(sprintf('h.pushbutton%d', pb(j))),'Enable','on')
end
pb = [11 12 13 21 22 23 31 32 33];
for j = 1:numel(pb)
    set(eval(sprintf('h.Q%d', pb(j))),'Enable','on')
end
set(h.full,'Enable', 'on');
set(h.edit50,'Enable', 'on');
set(h.edit50,'String', num2str(h.dat.cl.threshold));

% setup different views of GUI
h.dat.maxmap = 2;
ops = h.dat.ops;
if isfield(ops, 'mimg1') && ~isempty(ops.mimg1)
    h.dat.mimg(:,:,h.dat.maxmap) = ops.mimg1(ops.yrange, ops.xrange);
    h.dat.mimg_proc(:,:,h.dat.maxmap) = normalize_image(h.dat.mimg(:,:,h.dat.maxmap));
end
h.dat.mimg(:,:,5) = 0;

h.dat.maxmap = h.dat.maxmap + 1;
if isfield(ops, 'mimgRED') && ~isempty(ops.mimgRED)
    h.dat.mimg(:,:,h.dat.maxmap) = ops.mimgRED(ops.yrange, ops.xrange);
    h.dat.mimg_proc(:,:,h.dat.maxmap) = normalize_image(h.dat.mimg(:,:,h.dat.maxmap));
end
h.dat.maxmap = h.dat.maxmap + 1;
if isfield(ops, 'mimgREDcorrected') && ~isempty(ops.mimgREDcorrected)
    h.dat.mimg(:,:,h.dat.maxmap) = ops.mimgREDcorrected;
    h.dat.mimg_proc(:,:,h.dat.maxmap) = normalize_image(h.dat.mimg(:,:,h.dat.maxmap));
end
h.dat.maxmap = h.dat.maxmap + 1;
if isfield(ops, 'Vcorr') && ~isempty(ops.Vcorr)
    h.dat.mimg(:,:,h.dat.maxmap) = ops.Vcorr;
    h.dat.mimg_proc(:,:,h.dat.maxmap) = normalize_image(h.dat.mimg(:,:,h.dat.maxmap));
end

h.dat.procmap = 0;

h.dat.map = 1;
h.dat.F.trace = [];
for i = 1:length(h.dat.Fcell)
    h.dat.F.trace = cat(2, h.dat.F.trace, h.dat.Fcell{i});
end
if isfield(h.dat, 'FcellNeu')
    h.dat.F.neurop = [];
    for i = 1:length(h.dat.FcellNeu)
        h.dat.F.neurop = cat(2, h.dat.F.neurop, h.dat.FcellNeu{i});
    end    
    
else
   h.dat.F.neurop = zeros(size(h.dat.F.trace), 'single');
end
h.dat.plot_neu = 1;

redraw_fluorescence(h);
redraw_figure(h);

guidata(hObject,h)
end

function pushbutton2_Callback(hObject, eventdata, h)
% variance explained mask
V0 = reshape(h.dat.res.lambda, h.dat.cl.Ly, h.dat.cl.Lx);
V0 = .75 * V0/mean(V0(V0>1e-5));
h.dat.img0.V = V0;
h = buildLambdaValue(h);
guidata(hObject,h);
redraw_figure(h);

function pushbutton1_Callback(hObject, eventdata, h)
% unit vector mask
V0 = reshape(h.dat.res.lambda0, h.dat.cl.Ly, h.dat.cl.Lx);
h.dat.img0.V = .75 * V0 / mean(V0(V0>1e-6));
h = buildLambdaValue(h);
guidata(hObject,h);
redraw_figure(h);

function pushbutton20_Callback(hObject, eventdata, h)
% original default hue
h.dat.cl.rands   = h.dat.cl.rands_orig;
h.dat.img1.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
h.dat.img2.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
guidata(hObject,h);
redraw_figure(h);

function pushbutton84_Callback(hObject, eventdata, h)
% save proc file and rules file
h.dat.F.trace = [];
dat = h.dat;
save([h.dat.filename(1:end-4) '_proc.mat'], 'dat')
%
h.st0(:,1) = double([h.dat.stat.iscell]);
%
statLabels = [{'labels'} h.statLabels];
prior = h.prior;

st = cat(1, h.st, h.st0);save(h.dat.cl.fpath, 'st', 'statLabels', 'prior')


function figure1_ResizeFcn(hObject, eventdata, h)

function Q11_Callback(hObject, eventdata, h)
iy = 1; ix = 1;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

function Q12_Callback(hObject, eventdata, h)
iy = 1; ix = 2;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

function Q13_Callback(hObject, eventdata, h)
iy = 1; ix = 3;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

function Q21_Callback(hObject, eventdata, h)
iy = 2; ix = 1;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

function Q22_Callback(hObject, eventdata, h)
iy = 2; ix = 2;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

function Q23_Callback(hObject, eventdata, h)
iy = 2; ix = 3;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

function Q31_Callback(hObject, eventdata, h)
iy = 3; ix = 1;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

function Q32_Callback(hObject, eventdata, h)
iy = 3; ix = 2;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

function Q33_Callback(hObject, eventdata, h)
iy = 3; ix = 3;
quadrant(hObject, h, iy, ix);
paint_quadbutton(h, iy, ix);

function paint_quadbutton(h, iy, ix)
for j = 1:3
    for i = 1:3
        if h.quadvalue(j,i)==1
            set(h.(sprintf('Q%d%d', j,i)), 'BackgroundColor', [.4 .4 .4]); 
        end
    end
end
set(h.(sprintf('Q%d%d', iy,ix)), 'BackgroundColor', [.8 .0 .3]); 

function full_Callback(hObject, eventdata, h)
h.dat.ylim = [0 h.dat.cl.Ly];
h.dat.xlim = [0 h.dat.cl.Lx];
guidata(hObject,h);
if h.dat.map==1
    redraw_figure(h);
else
    redraw_meanimg(h);
end

function quadrant(hObject, h, iy, ix)
h.dat.ylim = [h.dat.figure.y0all(iy) h.dat.figure.y1all(iy+1)];
h.dat.xlim = [h.dat.figure.x0all(ix) h.dat.figure.x1all(ix+1)];
h.quadvalue(iy, ix) = 1;

guidata(hObject,h);

if h.dat.map==1
    redraw_figure(h);
else
    redraw_meanimg(h);
end

function figure1_WindowKeyPressFcn(hObject, eventdata, h)
switch eventdata.Key
    case 'f'
        % flip currently selected unit
                    h.dat.stat(h.dat.F.ichosen).iscell = 1 - ...
                h.dat.stat(h.dat.F.ichosen).iscell;
        
        h = splitROIleftright(h);
        h = buildLambdaValue(h);
        guidata(hObject,h);
        if h.dat.maxmap==1
            redraw_figure(h);
        end
    case 'q'
        pushbutton87_Callback(hObject, eventdata, h);
    case 'w'
        pushbutton89_Callback(hObject, eventdata, h);
    case 'e'
        if size(h.dat.mimg, 3)>2
            pushbutton90_Callback(hObject, eventdata, h);
        end
    case 'r'
        pushbutton92_Callback(hObject, eventdata, h);
    case 's'
        pushbutton103_Callback(hObject, eventdata, h);
    case 'p'
        pushbutton86_Callback(hObject, eventdata, h);
end

function figure1_WindowButtonDownFcn(hObject, eventdata, h)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
z = round(eventdata.Source.CurrentAxes.CurrentPoint(1,:));
x = round(z(1));
y  = round(z(2));

if x>=1 && y>=1 && x<=h.dat.cl.Lx && y<=h.dat.cl.Ly && h.dat.res.iclust(y,x)>0
    h.dat.F.ichosen = h.dat.res.iclust(y, x);
    ichosen = h.dat.F.ichosen;
    
    switch eventdata.Source.SelectionType
        case 'alt'
            % flip currently selected unit
            h.dat.stat(ichosen).iscell = 1 - ...
                h.dat.stat(ichosen).iscell;
            
            h = splitROIleftright(h);
            h = buildLambdaValue(h);
        case 'extend'
            h.dat.cl.redcell(ichosen) = 1 -  h.dat.cl.redcell(ichosen);
            
            if h.dat.cl.redcell(ichosen)==1
                h.dat.cl.rands(ichosen) = 0;
            else
                h.dat.cl.rands(ichosen) = h.dat.cl.rands_orig(ichosen);
            end
            h.dat.img1.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
            h.dat.img2.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
            
            if h.dat.cl.redcell(ichosen)
                display('red')
            else
                display('not red')
            end
    end
    
    redraw_fluorescence(h);
    
    Sat = ones(h.dat.cl.Ly, h.dat.cl.Lx);
    Sat(h.dat.res.iclust==h.dat.F.ichosen) = 0;
    h.dat.img1.Sat     = Sat;
    h.dat.img2.Sat     = Sat;
    h = buildLambdaValue(h);
    
    guidata(hObject,h);
    redraw_figure(h);
    
    str = [];
    for j =1:length(h.statLabels)
       if isfield(h.dat.stat, h.statLabels{j})
           sl = eval(sprintf('h.dat.stat(ichosen).%s', h.statLabels{j}));
           strnew = sprintf('%s = %2.2f \n', h.statLabels{j}, sl);
           str = cat(2, str, strnew);
       end
    end
    
   set(h.text54,'String', str);
end

function pushbutton86_Callback(hObject, eventdata, h)

h.dat.procmap = 1 -  h.dat.procmap;
if h.dat.map>1
    redraw_meanimg(h);
end

% if h.dat.map < h.dat.maxmap
%     h.dat.map = h.dat.map + 1;
%     redraw_meanimg(h);
% else
%     h.dat.map = 1;
%     redraw_figure(h);
% end
guidata(hObject,h);

function pushbutton87_Callback(hObject, eventdata, h)
 h.dat.map = 1;
redraw_figure(h);
guidata(hObject,h);

function pushbutton89_Callback(hObject, eventdata, h)
h.dat.map = 2;
redraw_meanimg(h);
guidata(hObject,h);

function pushbutton90_Callback(hObject, eventdata, h)
 h.dat.map = 3;
 redraw_meanimg(h);
 guidata(hObject,h);
 
% RED CORRECTED BUTTON
function pushbutton92_Callback(hObject, eventdata, h)
h.dat.map = 4;
redraw_meanimg(h);
guidata(hObject,h);

function pushbutton91_Callback(hObject, eventdata, h)
h.dat.plot_neu = 1 - h.dat.plot_neu;
redraw_fluorescence(h)
guidata(hObject,h);

function pushbutton93_Callback(hObject, eventdata, h)
rootS2p = which('run_pipeline');
if isempty(rootS2p)
    error('could not identify Suite2p location! where is new_main.m?')
end
rootS2p = fileparts(rootS2p);
rootS2p = fullfile(rootS2p, 'configFiles');

[filename1,filepath1]   = uigetfile(rootS2p, 'Select classifier file');
if filename1
    h.dat.cl.fpath          = fullfile(filepath1, filename1);
    h              = classROI(h);
    
    h = splitROIleftright(h);
    h = buildLambdaValue(h);
    redraw_figure(h);
    
%     hload = load(h.dat.cl.fpath);
%     h.st = hload.st;
%     h.save_cls = h.dat.cl.fpath;
%     h.prior = prior;

    guidata(hObject,h);
end

function edit50_Callback(hObject, eventdata, h)
h.dat.cl.threshold = str2double(get(h.edit50,'String'));
for j = 1:length(h.dat.stat)
    h.dat.stat(j).iscell = h.dat.stat(j).cellProb > h.dat.cl.threshold;
end
h = splitROIleftright(h);
h = buildLambdaValue(h);
redraw_figure(h);

guidata(hObject,h);

function edit50_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbutton96_Callback(hObject, eventdata, h)
hval = [h.dat.stat.skew];
h.dat.cl.rands   = .1 + .8 * min(1, hval/mean(hval));
h = buildHue(h);
h = buildLambdaValue(h);
redraw_figure(h);
guidata(hObject,h);

function pushbutton95_Callback(hObject, eventdata, h)
h.dat.cl.rands   = .1 + .8 * [h.dat.stat.cellProb];
h = buildHue(h);
h = buildLambdaValue(h);
redraw_figure(h);
guidata(hObject,h);

function pushbutton98_Callback(hObject, eventdata, h)
h.dat.cl.rands   = .1 + .8 * rand(1, h.dat.ops.Nk);
h = buildHue(h);
h = buildLambdaValue(h);
redraw_figure(h);
guidata(hObject,h);

function pushbutton99_Callback(hObject, eventdata, h)
hval = max(0, [h.dat.stat.cmpct]-1);
h.dat.cl.rands   = .1 + .8 * min(1, hval/mean(hval));
h = buildHue(h);
h = buildLambdaValue(h);
redraw_figure(h);
guidata(hObject,h);

function pushbutton100_Callback(hObject, eventdata, h)
hval = [h.dat.stat.footprint];
h.dat.cl.rands   = .1 + .8 * min(1, hval/mean(hval));
h = buildHue(h);
h = buildLambdaValue(h);
redraw_figure(h);
guidata(hObject,h);

function pushbutton101_Callback(hObject, eventdata, h)
% new classifier button
rootS2p = which('run_pipeline');
if isempty(rootS2p)
    warndlg('Suite2p location not in PATH! where is run_pipeline.m?')
    [filename1,rootS2p]=uigetfile(root, 'Select run_pipeline.m');
    
end
rootS2p = fileparts(rootS2p);
rootS2p = fullfile(rootS2p, 'configFiles');

def_name = fullfile(rootS2p, 'cl_new.mat');

[filename1,filepath1]   = uigetfile(fullfile(rootS2p, 'priors'), ...
        'Select priors file, will be used to seed the new classifier');
if filename1
    prior_file = fullfile(filepath1, filename1);
else
    error('you must choose a priors file')
end

[FileName,PathName] = uiputfile('*.mat', 'Create new classifier', def_name); 

if FileName
    load(prior_file);
    st = [];
    save(fullfile(PathName, FileName), 'st', 'prior', 'statLabels')
    
    h.dat.cl.fpath          = fullfile(PathName, FileName);
    h                       = classROI(h);
    
    h = splitROIleftright(h);
    h = buildLambdaValue(h);
    redraw_figure(h);
    
    hload = load(h.dat.cl.fpath);
    h.st = hload.st;
    h.save_cls = h.dat.cl.fpath;
    
    guidata(hObject,h);
else
    error('you did not create a new classifier file')
end

function pushbutton102_Callback(hObject, eventdata, h)
hval = [h.dat.stat.mimgProjAbs];
h.dat.cl.rands   = .1 + .8 * min(1, hval/mean(hval));
h = buildHue(h);
h = buildLambdaValue(h);
redraw_figure(h);
guidata(hObject,h);

function pushbutton103_Callback(hObject, eventdata, h)
 h.dat.map = 5;
redraw_meanimg(h);
guidata(hObject,h);

function pushbutton64_Callback(hObject, eventdata, h)

msgbox('Instructions for using this GUI! More detailed guide in included pdf \n You can keep this box open while interacting with the GUI. ', 'Instructions!');
