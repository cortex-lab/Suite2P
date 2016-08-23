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

h.output = hObject;
guidata(hObject, h);
% UIWAIT makes new_main wait for user response (see UIRESUME)
% uiwait(h.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = new_main_OutputFcn(hObject, eventdata, h) 
% varargout  cell array for returning output args (see VARARGOUT);
varargout{1} = h.output;

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
    
    h.dat.cl.Mrs      = [h.dat.stat.mrs]./[h.dat.stat.mrs0];
    h.dat.cl.npix     = [h.dat.stat.npix];
    h.dat.cl.Ly       = numel(h.dat.ops.yrange);
    h.dat.cl.Lx       = numel(h.dat.ops.xrange);
    h.dat.cl.MeanM    = 2*mean(h.dat.res.M);
    h.dat.cl.excluded_pixels  = zeros(h.dat.cl.Ly, h.dat.cl.Lx);
    h.dat.cl.excluded_regions = zeros(h.dat.cl.Ly, h.dat.cl.Lx);
    h.dat.cl.excl_pix_perc    = zeros(h.dat.cl.Ly, h.dat.cl.Lx);
    h.dat.cl.topregion        = ones(h.dat.cl.Ly, h.dat.cl.Lx);
%     if isfield(h.dat.stat, 'parent')
        h = get_parent_stats(h);
%     end
    
    h.dat.res.iclust = reshape(h.dat.res.iclust, h.dat.cl.Ly, h.dat.cl.Lx);
    
    Nk = h.dat.ops.Nk;
    h.dat.ops.Nk = numel(h.dat.stat);
    h.dat.cl.rands_orig   = .1 + .8 * rand(1, h.dat.ops.Nk);
    h.dat.cl.rands        = h.dat.cl.rands_orig;
    
    if isfield(h.dat, 'clustrules')
         % ROI rules
         h.dat.res.Mrs_thresh_orig = h.dat.clustrules.Compact;
         h.dat.cl.npix_low_orig    = h.dat.clustrules.MinNpix;
         h.dat.cl.npix_high_orig   = h.dat.clustrules.MaxNpix;
    else
        % ROI rules
        h.dat.res.Mrs_thresh_orig   = 3;
        h.dat.cl.npix_low_orig      = 20;
        h.dat.cl.npix_high_orig     = 500;
    end

    % parent rules
    h.dat.cl.mrs_parent_max = Inf;
    h.dat.cl.npix_res_max   = Inf;
    h.dat.cl.npix_par_max   = Inf;
    h.dat.cl.nreg_max       = Inf;
    h.dat.cl.VperPix_min    = 0;
    
    h = setOriginalThresh(h);
    
    set(h.edit35,'String', num2str(h.dat.res.Mrs_thresh));
    set(h.edit39,'String', num2str(h.dat.cl.npix_high));
    set(h.edit40,'String', num2str(h.dat.cl.npix_low));
    
    % set all quadrants as not visited
    h.quadvalue = zeros(3);
    for j = 1:3
        for i = 1:3
            set(h.(sprintf('Q%d%d', j,i)), 'BackgroundColor',[.92 .92 .92]);
        end
    end
    % start with unit vector map
    lam = h.dat.res.lambda;
    h.dat.img0.V = max(0, min(1, .5 * reshape(lam, h.dat.cl.Ly, h.dat.cl.Lx)/mean(lam(:))));
    
    h.dat.ylim = [0 h.dat.cl.Ly];
    h.dat.xlim = [0 h.dat.cl.Lx];    
    
    h.dat.cl.manual  = zeros(h.dat.ops.Nk, 1);
    if ~isfield(h.dat.cl,'redcell')
        h.dat.cl.redcell = zeros(h.dat.ops.Nk, 1);
    end
    h                = splitROIleftright(h);
    
    icell = find(h.dat.cl.iscell);
    if ~isempty(icell)
        h.dat.F.ichosen = icell(1); %ceil(rand * numel(icell))
    else
        h.dat.F.ichosen = 1; %ceil(rand * numel(icell))
    end
    
    Sat = ones(h.dat.cl.Ly, h.dat.cl.Lx);
    Sat(h.dat.res.iclust==h.dat.F.ichosen) = 0;
    h.dat.img1.Sat     = Sat;
    h.dat.img2.Sat     = Sat;
    
    h = buildHue(h);
    h = buildLambdaValue(h);
    
    % loop through redcells and set h.dat.cl.rands(h.dat.F.ichosen) = 0
    for j = find(h.dat.cl.redcell)
        h.dat.F.ichosen = j;
        h.dat.cl.rands(h.dat.F.ichosen) = 0;
    end
    if ~isempty(icell)
        h.dat.F.ichosen = icell(1); %ceil(rand * numel(icell))
    else
        h.dat.F.ichosen = 1; %ceil(rand * numel(icell))
    end
    h = buildHue(h);
    h = buildLambdaValue(h);
    
    % x and y limits on subquadrants
    h.dat.figure.x0all = round(linspace(0, 19/20*h.dat.cl.Lx, 4));
    h.dat.figure.y0all = round(linspace(0, 19/20*h.dat.cl.Ly, 4));
    h.dat.figure.x1all = round(linspace(1/20 * h.dat.cl.Lx, h.dat.cl.Lx, 4));
    h.dat.figure.y1all = round(linspace(1/20 * h.dat.cl.Ly, h.dat.cl.Ly, 4));
    
    h.dat.F.Fcell = h.dat.Fcell; h.dat.Fcell = [];    
    
    if isfield(h.dat, 'FcellNeu')
        h.dat.F.FcellNeu = h.dat.FcellNeu; h.dat.FcellNeu = [];
        if mean(sign(h.dat.F.FcellNeu{1}(:)))<0
            for j = 1:length(h.dat.F.FcellNeu)
                h.dat.F.FcellNeu{j} = - h.dat.F.FcellNeu{j};
                h.dat.F.Fcell{j} = h.dat.F.Fcell{j} + h.dat.F.FcellNeu{j};
            end
        end    
        
%         if isfield(h.dat.cl, 'dcell')
%             for k = 1:length(h.dat.cl.dcell)
%                 for j = 1:length(h.dat.F.FcellNeu)
%                     if isfield(h.dat.cl.dcell{k}, 'B')
%                         c2 = h.dat.cl.dcell{k}.B(3);
%                         c1 = h.dat.cl.dcell{k}.B(2);
%                         h.dat.F.FcellNeu{j}(k+Nk, :) = c1 + c2 * h.dat.F.FcellNeu{j}(k+Nk, :);
%                     end
%                 end
%             end
%         end
    end
end

h.dat.maxmap = 1;
ops = h.dat.ops;
if isfield(ops, 'mimg1') && ~isempty(ops.mimg1)
    h.dat.maxmap = h.dat.maxmap + 1;
    h.dat.mimg(:,:,h.dat.maxmap) = ops.mimg1(ops.yrange, ops.xrange);
    h.dat.mimg_proc(:,:,h.dat.maxmap) = normalize_image(h.dat.mimg(:,:,h.dat.maxmap));
end
if isfield(ops, 'mimgRED') && ~isempty(ops.mimgRED)
    h.dat.maxmap = h.dat.maxmap + 1;
    h.dat.mimg(:,:,h.dat.maxmap) = ops.mimgRED(ops.yrange, ops.xrange);
    h.dat.mimg_proc(:,:,h.dat.maxmap) = normalize_image(h.dat.mimg(:,:,h.dat.maxmap));
end
if isfield(ops, 'mimgREDcorrected') && ~isempty(ops.mimgREDcorrected)
    h.dat.maxmap = h.dat.maxmap + 1;
    h.dat.mimg(:,:,h.dat.maxmap) = ops.mimgREDcorrected;
    h.dat.mimg_proc(:,:,h.dat.maxmap) = normalize_image(h.dat.mimg(:,:,h.dat.maxmap));
end

h.dat.procmap = 0;

h.dat.map = 1;
h.dat.F.trace = [];
for i = 1:length(h.dat.F.Fcell)
    h.dat.F.trace = cat(2, h.dat.F.trace, h.dat.F.Fcell{i});
end
if isfield(h.dat.F, 'FcellNeu')
    h.dat.F.neurop = [];
    for i = 1:length(h.dat.F.FcellNeu)
        h.dat.F.neurop = cat(2, h.dat.F.neurop, h.dat.F.FcellNeu{i});
    end    
    
else
   h.dat.F.neurop = zeros(size(h.dat.F.trace), 'single');
end
h.dat.plot_neu = 0;

redraw_fluorescence(h);
redraw_figure(h);

guidata(hObject,h)
end

function pushbutton61_Callback(hObject, eventdata, h)
% keep TOP variance region
h.dat.cl.topregion = zeros(h.dat.cl.Ly, h.dat.cl.Lx);
xs = repmat(1:h.dat.cl.Lx, h.dat.cl.Ly, 1);
ys = repmat((1:h.dat.cl.Ly)', 1, h.dat.cl.Lx);

for k = 1:length(h.dat.stat)
    if ~isempty(h.dat.stat(k).Vregion)
        [~, itop] = max(h.dat.stat(k).Vregion);
        h.dat.cl.topregion(h.dat.stat(k).region{itop}) = 1;
        h.dat.cl.npix(k) = h.dat.stat(k).npixels(itop);
    end
end

h = buildLambdaValue(h);
redraw_figure(h);
guidata(hObject,h)

function pushbutton18_Callback(hObject, eventdata, h)
function slider5_Callback(hObject, eventdata, h)
function slider5_CreateFcn(hObject, eventdata, h)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function slider4_Callback(hObject, eventdata, h)
function slider4_CreateFcn(hObject, eventdata, h)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function slider2_Callback(hObject, eventdata, h)
function slider2_CreateFcn(hObject, eventdata, h)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
    
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
function slider1_CreateFcn(hObject, eventdata, h)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function pushbutton3_Callback(hObject, eventdata, h)
% binary mask
h.dat.img0.V = ones(h.dat.cl.Ly, h.dat.cl.Lx);
h = buildLambdaValue(h);

guidata(hObject,h);
redraw_figure(h);

function pushbutton2_Callback(hObject, eventdata, h)
% variance explained mask
h.dat.img0.V = reshape(h.dat.res.M, h.dat.cl.Ly, h.dat.cl.Lx)/h.dat.cl.MeanM;
h = buildLambdaValue(h);
guidata(hObject,h);
redraw_figure(h);

function pushbutton1_Callback(hObject, eventdata, h)
% unit vector mask
h.dat.img0.V = 10 * reshape(h.dat.res.lambda, h.dat.cl.Ly, h.dat.cl.Lx);
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

function pushbutton4_Callback(hObject, eventdata, h)
% randomize hue
rng('shuffle') 
h.dat.cl.rands     = rand(1, h.dat.ops.Nk);
h.dat.cl.rands(1)  = .15;
h.dat.img1.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
h.dat.img2.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);

guidata(hObject,h);
redraw_figure(h);

function pushbutton36_Callback(hObject, eventdata, h)
function pushbutton34_Callback(hObject, eventdata, h)
function popupmenu4_Callback(hObject, eventdata, h)
function popupmenu4_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu5_Callback(hObject, eventdata, h)
function popupmenu5_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu6_Callback(hObject, eventdata, h)
function popupmenu6_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu7_Callback(hObject, eventdata, h)
function popupmenu7_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function slider8_Callback(hObject, eventdata, h)
function slider8_CreateFcn(hObject, eventdata, h)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function pushbutton35_Callback(hObject, eventdata, h)

function pushbutton37_Callback(hObject, eventdata, h)
function pushbutton43_Callback(hObject, eventdata, h)
function pushbutton43_ButtonDownFcn(hObject, eventdata, h)
function pushbutton44_Callback(hObject, eventdata, h)
function pushbutton45_Callback(hObject, eventdata, h)
function pushbutton46_Callback(hObject, eventdata, h)
function pushbutton47_Callback(hObject, eventdata, h)
function edit2_Callback(hObject, eventdata, h)
function edit2_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton38_Callback(hObject, eventdata, h)
function pushbutton39_Callback(hObject, eventdata, h)
function pushbutton40_Callback(hObject, eventdata, h)
function pushbutton41_Callback(hObject, eventdata, h)
function pushbutton42_Callback(hObject, eventdata, h)
function edit1_Callback(hObject, eventdata, h)
function edit1_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton57_Callback(hObject, eventdata, h)
function edit3_Callback(hObject, eventdata, h)
function edit3_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit4_Callback(hObject, eventdata, h)
function edit4_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit5_Callback(hObject, eventdata, h)
function edit5_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit6_Callback(hObject, eventdata, h)
function edit6_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton50_Callback(hObject, eventdata, h)
function pushbutton48_Callback(hObject, eventdata, h)
function pushbutton54_Callback(hObject, eventdata, h)
function pushbutton55_Callback(hObject, eventdata, h)
function pushbutton56_Callback(hObject, eventdata, h)
function pushbutton52_Callback(hObject, eventdata, h)
function pushbutton53_Callback(hObject, eventdata, h)
function edit7_Callback(hObject, eventdata, h)
function edit7_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit8_Callback(hObject, eventdata, h)
function edit8_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit10_Callback(hObject, eventdata, h)
function edit10_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit11_Callback(hObject, eventdata, h)
function edit11_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton59_Callback(hObject, eventdata, h)
function edit12_Callback(hObject, eventdata, h)
function edit12_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit29_Callback(hObject, eventdata, h)
function edit29_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit30_Callback(hObject, eventdata, h)
function edit30_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit31_Callback(hObject, eventdata, h)
function edit31_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit32_Callback(hObject, eventdata, h)
function edit32_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbutton63_Callback(hObject, eventdata, h)
function edit33_Callback(hObject, eventdata, h)
h.dat.cl.pixthresh_percent = str2double(get(h.edit33,'String'));
h = exclude_pixels_percent(h);
h = buildLambdaValue(h);
redraw_figure(h);
guidata(hObject,h);

function h = exclude_pixels_percent(h)
h.dat.cl.excl_pix_perc = zeros(h.dat.cl.Ly, h.dat.cl.Lx);
for k = 1:h.dat.ops.Nk
    which_pix = find(h.dat.res.iclust==k);
   [Msort, isort] = sort(h.dat.res.M(which_pix), 'ascend'); 
   Msort = cumsum(Msort);
   Msort = 100 * Msort/max(Msort);
   ifi = find(Msort>h.dat.cl.pixthresh_percent, 1);
   h.dat.cl.excl_pix_perc(which_pix(isort(1:ifi))) = 1;
end

function edit33_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit34_Callback(hObject, eventdata, h)
function edit34_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton67_Callback(hObject, eventdata, h)
function edit37_Callback(hObject, eventdata, h)
function edit37_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton70_Callback(hObject, eventdata, h)
function pushbutton65_Callback(hObject, eventdata, h)

function h = pushbutton64_Callback(hObject, eventdata, h)
[x, y] = ginput(1);
x = min(max(1, round(x)), h.dat.cl.Lx); 
y = min(max(1, round(y)), h.dat.cl.Ly); 
h.dat.F.ichosen = h.dat.res.iclust(y, x);
redraw_fluorescence(h);

Sat = ones(h.dat.cl.Ly, h.dat.cl.Lx);
Sat(h.dat.res.iclust==h.dat.F.ichosen) = 0;
h.dat.img1.Sat     = Sat;
h.dat.img2.Sat     = Sat;
h = buildLambdaValue(h);

guidata(hObject,h);
redraw_figure(h);

function edit36_Callback(hObject, eventdata, h)
function edit36_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit35_Callback(hObject, eventdata, h)
h.dat.res.Mrs_thresh = str2double(get(h.edit35,'String'));
h = splitROIleftright(h);
h = buildLambdaValue(h);
redraw_figure(h);
guidata(hObject,h);

function edit35_CreateFcn(hObject, eventdata, h)
set(hObject,'String', num2str(2));
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton76_Callback(hObject, eventdata, h)
function pushbutton77_Callback(hObject, eventdata, h)
function pushbutton78_Callback(hObject, eventdata, h)
function pushbutton73_Callback(hObject, eventdata, h)
function pushbutton74_Callback(hObject, eventdata, h)
function pushbutton75_Callback(hObject, eventdata, h)
function pushbutton83_Callback(hObject, eventdata, h)
function pushbutton84_Callback(hObject, eventdata, h)
h.dat.F.trace = [];
dat = h.dat;
save([h.dat.filename(1:end-4) '_proc.mat'], 'dat')

function pushbutton79_Callback(hObject, eventdata, h)
function pushbutton80_Callback(hObject, eventdata, h)

function edit39_Callback(hObject, eventdata, h)
h.dat.cl.npix_high = str2double(get(h.edit39,'String'));
h = splitROIleftright(h);
h = buildLambdaValue(h);
redraw_figure(h);
guidata(hObject,h);

function edit39_CreateFcn(hObject, eventdata, h)
set(hObject,'String', num2str(400));
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit40_Callback(hObject, eventdata, h)
h.dat.cl.npix_low = str2double(get(h.edit40,'String'));
h = splitROIleftright(h);
h = buildLambdaValue(h);
redraw_figure(h);
guidata(hObject,h);

function edit40_CreateFcn(hObject, eventdata, h)
set(hObject,'String', num2str(30));
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu11_Callback(hObject, eventdata, h)
function popupmenu11_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit41_Callback(hObject, eventdata, h)
function edit41_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit42_Callback(hObject, eventdata, h)
function edit42_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function figure1_ResizeFcn(hObject, eventdata, h)

function edit43_Callback(hObject, eventdata, h)
h.dat.cl.pixthresh_var = str2double(get(h.edit43,'String'));
h = excluded_pixels(h);
h = buildLambdaValue(h);
redraw_figure(h);
guidata(hObject,h);


function edit43_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit44_Callback(hObject, eventdata, h)
function edit44_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in Q11.
function Q11_Callback(hObject, eventdata, h)
iy = 1; ix = 1;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

% --- Executes on button press in Q12.
function Q12_Callback(hObject, eventdata, h)
iy = 1; ix = 2;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

% --- Executes on button press in Q13.
function Q13_Callback(hObject, eventdata, h)
iy = 1; ix = 3;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

% --- Executes on button press in Q21.
function Q21_Callback(hObject, eventdata, h)
iy = 2; ix = 1;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

% --- Executes on button press in Q22.
function Q22_Callback(hObject, eventdata, h)
iy = 2; ix = 2;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

% --- Executes on button press in Q23.
function Q23_Callback(hObject, eventdata, h)
iy = 2; ix = 3;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

% --- Executes on button press in Q31.
function Q31_Callback(hObject, eventdata, h)
iy = 3; ix = 1;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

% --- Executes on button press in Q32.
function Q32_Callback(hObject, eventdata, h)
iy = 3; ix = 2;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

% --- Executes on button press in Q33.
function Q33_Callback(hObject, eventdata, h)
iy = 3; ix = 3;
quadrant(hObject, h, iy, ix);
paint_quadbutton(h, iy, ix);

function paint_quadbutton(h, iy, ix)
for j = 1:3
    for i = 1:3
        if h.quadvalue(j,i)==1
            set(h.(sprintf('Q%d%d', j,i)), 'BackgroundColor','yellow'); 
        end
    end
end
set(h.(sprintf('Q%d%d', iy,ix)), 'BackgroundColor','red'); 

% --- Executes on button press in full.
function full_Callback(hObject, eventdata, h)
h.dat.ylim = [0 h.dat.cl.Ly];
h.dat.xlim = [0 h.dat.cl.Lx];
guidata(hObject,h);
redraw_figure(h);

function quadrant(hObject, h, iy, ix)
h.dat.ylim = [h.dat.figure.y0all(iy) h.dat.figure.y1all(iy+1)];
h.dat.xlim = [h.dat.figure.x0all(ix) h.dat.figure.x1all(ix+1)];
h.quadvalue(iy, ix) = 1;

guidata(hObject,h);
redraw_figure(h);

% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, h)
switch eventdata.Key
    case 'f'
        % flip currently selected unit
        h.dat.cl.manual(h.dat.F.ichosen) = .5 - h.dat.cl.iscell(h.dat.F.ichosen);
        h = splitROIleftright(h);
        h = buildLambdaValue(h);
        guidata(hObject,h);
        if h.dat.maxmap==1
            redraw_figure(h);
        end
    case 's'
        % manual selection of units
        pushbutton64_Callback(hObject, eventdata, h);
    case 'q'
        h.dat.map = 1;
        pushbutton87_Callback(hObject, eventdata, h);
    case 'w'
        h.dat.map = 2;
        pushbutton89_Callback(hObject, eventdata, h);
    case 'e'
        h.dat.map = 3;
        if h.dat.maxmap>2
            pushbutton90_Callback(hObject, eventdata, h);
        end
    case 'r'
        h.dat.map = 3;
        if h.dat.maxmap>3
            pushbutton92_Callback(hObject, eventdata, h);
        end
        
    case 'p'
        pushbutton86_Callback(hObject, eventdata, h);
end


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, h)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
z = round(eventdata.Source.CurrentAxes.CurrentPoint(1,:));
x = round(z(1));
y  = round(z(2));

% disp(eventdata.Source.SelectionType)
% keyboard;
% manual selection of units
x = min(max(1, round(x)), h.dat.cl.Lx);
y = min(max(1, round(y)), h.dat.cl.Ly);
h.dat.F.ichosen = h.dat.res.iclust(y, x);

switch eventdata.Source.SelectionType
    case 'alt'
        % flip currently selected unit
        h.dat.cl.manual(h.dat.F.ichosen) = .5 - h.dat.cl.iscell(h.dat.F.ichosen);
        h = splitROIleftright(h);
        h = buildLambdaValue(h);
    case 'open'
        % unpin the manual selection on this cell
        h.dat.cl.manual(h.dat.F.ichosen) = 0;
        h = splitROIleftright(h);
        h = buildLambdaValue(h);
    case 'extend'
         h.dat.cl.redcell(h.dat.F.ichosen) = 1 -  h.dat.cl.redcell(h.dat.F.ichosen);
        
        if h.dat.cl.redcell(h.dat.F.ichosen)==1
            h.dat.cl.rands(h.dat.F.ichosen) = 0;
        else
             h.dat.cl.rands(h.dat.F.ichosen) = h.dat.cl.rands_orig(h.dat.F.ichosen);
        end
        h.dat.img1.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
        h.dat.img2.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
        
        if h.dat.cl.redcell(h.dat.F.ichosen)
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



function edit45_Callback(hObject, eventdata, h)
h.dat.cl.mrs_parent_max = str2double(get(h.edit45,'String'));
h = splitROIleftright(h);
h = buildLambdaValue(h);
redraw_figure(h);
guidata(hObject,h);

% --- Executes during object creation, after setting all properties.
function edit45_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject,'String', num2str(Inf));
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit46_Callback(hObject, eventdata, h)
h.dat.cl.npix_par_max = str2double(get(h.edit46,'String'));
h = splitROIleftright(h);
h = buildLambdaValue(h);
redraw_figure(h);
guidata(hObject,h);

% --- Executes during object creation, after setting all properties.
function edit46_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'String', num2str(Inf));
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit47_Callback(hObject, eventdata, h)
h.dat.cl.npix_res_max = str2double(get(h.edit47,'String'));
h = splitROIleftright(h);
h = buildLambdaValue(h);
redraw_figure(h);
guidata(hObject,h);


% --- Executes during object creation, after setting all properties.
function edit47_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject,'String', num2str(Inf));
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit48_Callback(hObject, eventdata, h)
h.dat.cl.nreg_max = str2double(get(h.edit48,'String'));
h = splitROIleftright(h);
h = buildLambdaValue(h);
redraw_figure(h);
guidata(hObject,h);

% --- Executes during object creation, after setting all properties.
function edit48_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject,'String', num2str(Inf));
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit49_Callback(hObject, eventdata, h)
h.dat.cl.VperPix_min = str2double(get(h.edit49,'String'));
h = splitROIleftright(h);
h = buildLambdaValue(h);
redraw_figure(h);
guidata(hObject,h);



% --- Executes during object creation, after setting all properties.
function edit49_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'String', num2str(0));
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton86.
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


% --- Executes on button press in pushbutton87.
function pushbutton87_Callback(hObject, eventdata, h)
 h.dat.map = 1;
redraw_figure(h);
guidata(hObject,h);


% --- Executes on button press in pushbutton89.
function pushbutton89_Callback(hObject, eventdata, h)
h.dat.map = 2;
redraw_meanimg(h);
guidata(hObject,h);

% --- Executes on button press in pushbutton90.
function pushbutton90_Callback(hObject, eventdata, h)
 h.dat.map = 3;
redraw_meanimg(h);
guidata(hObject,h);

% RED CORRECTED BUTTON
function pushbutton92_Callback(hObject, eventdata, h)
h.dat.map = 4;
redraw_meanimg(h);
guidata(hObject,h);






% --- Executes on button press in pushbutton91.
function pushbutton91_Callback(hObject, eventdata, h)
h.dat.plot_neu = 1 - h.dat.plot_neu;
redraw_fluorescence(h)
guidata(hObject,h);



% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton86.
function pushbutton86_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton86 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

