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
axes(h.axes5);
set(gca, 'xtick', [], 'ytick', [])
%set(gca, 'xcolor', 'w', 'ycolor', 'w')
ylabel('scale','fontweight','bold','fontsize',8);

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
    [filename1,filepath1]=uigetfile(fullfile(root, 'F*.mat'), 'Select Data File');
    set(h.figure1, 'Name', filename1);
    
    % construct dat. Everything is loaded except F and Fneu.
    h.dat = load(fullfile(filepath1, filename1));
    
    flag = 1;
catch
end

if flag
    % if the user selected a file, do all the initializations
rng('default')

% keyboard;
init = 0;
if isfield(h.dat, 'dat')
    h.dat = h.dat.dat;
    if isfield(h.dat, 'cl')
      h                   = identify_classifier(h);    
      h                   = classROI(h);
      init = 1;
    else
      init = 0;
    end
end

if init==0
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
%     h.dat.res.lambda = reshape(h.dat.res.lambda, h.dat.cl.Ly, h.dat.cl.Lx);
    
    h.dat.ops.Nk = numel(h.dat.stat);
    h.dat.cl.rands_orig   = .1 + .65 * rand(1, h.dat.ops.Nk);
    h.dat.cl.rands        = h.dat.cl.rands_orig;
    
    if isfield(h.dat.ops, 'clustrules')
       h.dat.clustrules = h.dat.ops.clustrules; 
    end
    
    % set up classifier
    h.dat.cl.threshold  = 0.5;
    h                   = identify_classifier(h);    
    h                   = classROI(h);
    
    % set all quadrants as not visited
    h.quadvalue = zeros(3);
    for j = 1:3
        for i = 1:3
            set(h.(sprintf('Q%d%d', j,i)), 'BackgroundColor',[.92 .92 .92]);
        end
    end
    
    h.dat.ylim = [0 h.dat.cl.Ly];
    h.dat.xlim = [0 h.dat.cl.Lx];    
    
    if ~isfield(h.dat.stat,'redcell')
        for j = 1:numel(h.dat.stat)
            h.dat.stat(j).redcell = 0;
            h.dat.stat(j).redprob = 0;
        end
    end    
   
    h.dat.F.ichosen = 1;
    
    % loop through redcells and set h.dat.cl.rands(h.dat.F.ichosen) = 0
    for j = find([h.dat.stat.redcell])
        h.dat.cl.rands(j) = 0;
    end
   
    % x and y limits on subquadrants
    h.dat.figure.x0all = round(linspace(0, 19/20*h.dat.cl.Lx, 4));
    h.dat.figure.y0all = round(linspace(0, 19/20*h.dat.cl.Ly, 4));
    h.dat.figure.x1all = round(linspace(1/20 * h.dat.cl.Lx, h.dat.cl.Lx, 4));
    h.dat.figure.y1all = round(linspace(1/20 * h.dat.cl.Ly, h.dat.cl.Ly, 4));
    
end

% activate all pushbuttons
pb = [84 93 101 86 87 89 90 92 103 98 95 96 102 99 100 1 2 104 112];
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
set(h.edit52,'Enable', 'on');
set(h.edit54,'Enable', 'on');
set(h.edit52,'String','-Inf');
set(h.edit54,'String','Inf');
h.statstr = 'npix';
h.statstrs = {'npix','cmpct','aspect_ratio','skew','std','footprint','mimgProj'};
h.statnum = 1;
h.statmins = -Inf*ones(1,7);
h.statmaxs = Inf*ones(1,7);
set_Bcolor(h, 1);
set_maskCcolor(h, 1);
% select unit normalized ROI brightness
h.dat.cl.vmap = 'unit';
set_maskBcolor(h, 1);
set(h.full, 'BackgroundColor', [1 0 0])

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
elseif isfield(ops, 'AlignToRedChannel') && ops.AlignToRedChannel == 1 && ...
        isfield(ops, 'mimg') && ~isempty(ops.mimg)
    h.dat.mimg(:,:,h.dat.maxmap) = ops.mimg(ops.yrange, ops.xrange);
    h.dat.mimg_proc(:,:,h.dat.maxmap) = normalize_image(h.dat.mimg(:,:,h.dat.maxmap));
end
h.dat.maxmap = h.dat.maxmap + 1;
if isfield(ops, 'mimgREDcorrected') && ~isempty(ops.mimgREDcorrected)
    if sum(size(ops.mimgREDcorrected)==[ops.Ly ops.Lx]) == 2
        h.dat.mimg(:,:,h.dat.maxmap) = ops.mimgREDcorrected(ops.yrange, ops.xrange);
    else
        h.dat.mimg(:,:,h.dat.maxmap) = ops.mimgREDcorrected;
    end
    h.dat.mimg_proc(:,:,h.dat.maxmap) = normalize_image(h.dat.mimg(:,:,h.dat.maxmap));
end
h.dat.maxmap = h.dat.maxmap + 1;
if isfield(ops, 'Vcorr') && ~isempty(ops.Vcorr)
    h.dat.mimg(:,:,h.dat.maxmap) = ops.Vcorr;
    h.dat.mimg_proc(:,:,h.dat.maxmap) = normalize_image(h.dat.mimg(:,:,h.dat.maxmap));
end

h.dat.procmap = 0;
h.dat.map = 1;

redraw_fluorescence(h);
redraw_figure(h);

guidata(hObject,h)
end

% --------------- MASK BRIGHTNESS ----------------%
function set_maskBcolor(h, ih)
set(h.pushbutton1, 'BackgroundColor', .94 * [1 1 1]); 
set(h.pushbutton2, 'BackgroundColor', .94 * [1 1 1]); 

switch ih
    case 1
        set(h.pushbutton2, 'BackgroundColor', [1 0 0]); 
    case 2
        set(h.pushbutton1, 'BackgroundColor', [1 0 0]); 
end

function pushbutton2_Callback(hObject, eventdata, h)
% variance explained mask
h.dat.cl.vmap = 'unit';
set_maskBcolor(h, 1)

redraw_figure(h);
guidata(hObject,h);

function pushbutton1_Callback(hObject, eventdata, h)
% unit vector mask
h.dat.cl.vmap = 'var';
set_maskBcolor(h, 2)

redraw_figure(h);
guidata(hObject,h);

% --- save proc file
function pushbutton84_Callback(hObject, eventdata, h)
% save proc file and rules file
h.dat.F.trace = [];
h=classifierFig(h);
h=skewFig(h);
h=meanimgFig(h);
h=cmpctFig(h);
h=footprintFig(h);
h=redFig(h);
h=ellipseFig(h);

dat = h.dat;
save([h.dat.filename(1:end-4) '_proc.mat'], 'dat')
%
h.st0(:,1) = double([h.dat.stat.iscell]);
%
statLabels  = h.statLabels;
prior       = h.prior;
st          = cat(1, h.st, h.st0);
save(h.dat.cl.fpath, 'st', 'statLabels', 'prior')


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
set(h.full, 'BackgroundColor', .92 * [1 1 1])

for j = 1:3
    for i = 1:3
        if h.quadvalue(j,i)==1
            set(h.(sprintf('Q%d%d', j,i)), 'BackgroundColor', [.4 .4 .4]); 
        end
    end
end
set(h.(sprintf('Q%d%d', iy,ix)), 'BackgroundColor', [1 0 0]); 

function full_Callback(hObject, eventdata, h)
h.dat.ylim = [0 h.dat.cl.Ly];
h.dat.xlim = [0 h.dat.cl.Lx];
if h.dat.map==1
    redraw_figure(h);
else
    redraw_meanimg(h);
end
set(h.full, 'BackgroundColor', [1 0 0]);
for i = 1:3
    for j = 1:3
        if h.(sprintf('Q%d%d', j,i)).BackgroundColor(1) >.99 
            set(h.(sprintf('Q%d%d', j,i)), 'BackgroundColor', [.4 .4 .4])
        end
    end
end
guidata(hObject,h);

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
        if h.dat.maxmap==1
            redraw_figure(h);
        end
        guidata(hObject,h);
    case 'q'
        pushbutton87_Callback(hObject, eventdata, h);
    case 'e'
        pushbutton89_Callback(hObject, eventdata, h);
    case 'r'
        if size(h.dat.mimg, 3)>2
            pushbutton90_Callback(hObject, eventdata, h);
        end
    case 't'
        pushbutton92_Callback(hObject, eventdata, h);
    case 'w'
        pushbutton103_Callback(hObject, eventdata, h);
    case 'p'
        pushbutton86_Callback(hObject, eventdata, h);
    case 'a'
        pushbutton98_Callback(hObject, eventdata, h);
    case 's'
        pushbutton95_Callback(hObject, eventdata, h);
    case 'd'
        pushbutton96_Callback(hObject, eventdata, h);
    case 'g'
        pushbutton112_Callback(hObject, eventdata, h);
    case 'z'
        pushbutton102_Callback(hObject, eventdata, h);
    case 'x'
        pushbutton99_Callback(hObject, eventdata, h);
    case 'c'
        pushbutton100_Callback(hObject, eventdata, h);
    case 'v'
        pushbutton104_Callback(hObject, eventdata, h);
end

% ------------------ CELL CLICKING!! -------------------------%
function figure1_WindowButtonDownFcn(hObject, eventdata, h)
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
        case 'extend'
            h.dat.stat(ichosen).redcell = 1 -  h.dat.stat(ichosen).redcell;
            
            if h.dat.stat(ichosen).redcell ==1
                h.dat.cl.rands(ichosen) = 0;
            else
                h.dat.cl.rands(ichosen) = h.dat.cl.rands_orig(ichosen);
            end
%             h.dat.img1.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
%             h.dat.img2.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
            
            if h.dat.stat(ichosen).redcell
                display('red')
            else
                display('not red')
            end
    end
    
    redraw_fluorescence(h);
    
%     Sat = ones(h.dat.cl.Ly, h.dat.cl.Lx);
%     Sat(h.dat.res.iclust==h.dat.F.ichosen) = 0;
%     h.dat.img1.Sat     = Sat;
%     h.dat.img2.Sat     = Sat;
%     h = buildLambdaValue(h);
    
    redraw_figure(h);
    guidata(hObject,h);
    
    str = sprintf('cell #%d\n',ichosen);
    
    labels = [h.statLabels(2:end), {'iscell'}, {'redcell'},{'redprob'}];
    for j =1:length(labels)
       if isfield(h.dat.stat, labels{j})
           sl = eval(sprintf('h.dat.stat(ichosen).%s', labels{j}));
           strnew = sprintf('%s = %2.2f \n', labels{j}, sl);
           str = cat(2, str, strnew);
       end
    end
    
   set(h.text54,'String', str);
end




%------------------------ BACKGROUND --------------------------------%

function set_Bcolor(h, ih)
pb = [87 103 89 90 92];
 
for j = 1:length(pb)
    if j==ih
        set(h.(sprintf('pushbutton%d', pb(ih))), 'BackgroundColor', [1 0 0]); 
    else
        if h.(sprintf('pushbutton%d', pb(j))).BackgroundColor(1)>.99
            set(h.(sprintf('pushbutton%d', pb(j))), 'BackgroundColor', .94 * [1 1 1]);
        end
    end
end

% --- roi
function pushbutton87_Callback(hObject, eventdata, h)
 h.dat.map = 1;
redraw_figure(h);
set_Bcolor(h, 1);
guidata(hObject,h);

% --- mean img
function pushbutton89_Callback(hObject, eventdata, h)
h.dat.map = 2;
redraw_meanimg(h);
set_Bcolor(h, 3);
guidata(hObject,h);

% --- mean red img
function pushbutton90_Callback(hObject, eventdata, h)
 h.dat.map = 3;
 redraw_meanimg(h);
set_Bcolor(h, 4);
 guidata(hObject,h);
 
% --- mean red - green img
function pushbutton92_Callback(hObject, eventdata, h)
h.dat.map = 4;
redraw_meanimg(h);
set_Bcolor(h, 5);
guidata(hObject,h);

% --- correlation map
function pushbutton103_Callback(hObject, eventdata, h)
h.dat.map = 5;
redraw_meanimg(h);
set_Bcolor(h, 2);
guidata(hObject,h);

% --- proc - processed images
function pushbutton86_Callback(hObject, eventdata, h)
h.dat.procmap = 1 -  h.dat.procmap;
if h.dat.map>1
    redraw_meanimg(h);
end

if h.dat.procmap>0
    set(h.pushbutton86, 'BackgroundColor', [1 0 0]); 
else
    set(h.pushbutton86, 'BackgroundColor', .94 * [1 1 1]); 
end
guidata(hObject,h);


%--------------------- MASK COLORS ----------------------%

function set_maskCcolor(h, ih)
pb = [98 95 96 102 99 100 104 112]; 
 
set_Bcolor(h, 1)

% create colormap for mask colors
axes(h.axes5);
cla;
hold off;
% if not random masks

if ih > 1
    cmap = h.dat.cl.cmap;
    cmap = cmap(~isnan(cmap(:,1)),:);
    [cmax,im] = max(cmap(:,1));
    hmax = cmap(im,2);
    [cmin,im] = min(cmap(:,1));
    hmin = cmap(im,2);

    nticks = 4;
    hrange = linspace(hmin, hmax, 100);
    hticks = hrange(round(linspace(1,length(hrange),nticks)));
    for j = 1:nticks
        [~, itick] = min(abs(cmap(:,2) - hticks(j)));
        xticklabels{j} = sprintf('%2.2f', cmap(itick,1));
    end
    I = hsv2rgb(cat(3, hrange', ones(100,1), ones(100,1)));
    imagesc(permute(I,[2 1 3]));
    set(gca,'ytick','');
    set(gca,'xtick', round(linspace(1,100,nticks)), 'xticklabel', xticklabels, ...
        'fontsize', 8, 'fontweight','bold')
    ylabel('scale','fontweight','bold','fontsize',8);
else
    set(gca,'ytick','','xtick','');
    ylabel('scale','fontweight','bold','fontsize',8);
end
    
for j = 1:length(pb)
    if j==ih
        set(h.(sprintf('pushbutton%d', pb(ih))), 'BackgroundColor', [1 0 0]); 
    else
        if h.(sprintf('pushbutton%d', pb(j))).BackgroundColor(1)>.99
            set(h.(sprintf('pushbutton%d', pb(j))), 'BackgroundColor', .94 * [1 1 1]);
        end
    end
end

% --- random
function pushbutton98_Callback(hObject, eventdata, h)
if isfield(h.dat.cl, 'rands_orig')
    h.dat.cl.rands = h.dat.cl.rands_orig;
else
    h.dat.cl.rands   = .1 + .7*rand(length(h.dat.stat), 1);
end
h.dat.cl.rands(logical([h.dat.stat.redcell])) = 0;
I = redraw_figure(h);
set_maskCcolor(h, 1);
guidata(hObject,h);

% --- classifier
function pushbutton95_Callback(hObject, eventdata, h)
h=classifierFig(h);
guidata(hObject,h);

function h = classifierFig(h)
hval0            = [h.dat.stat.cellProb];
hval             = .6 * (1 - hval0) + .15;
h.dat.cl.rands   = hval;
h.dat.cl.cmap    = [hval0(:) hval(:)]; 
I = redraw_figure(h);
h.dat.cl.classifierFig = I;
set_maskCcolor(h, 2);

% --- skew
function pushbutton96_Callback(hObject, eventdata, h)
h = skewFig(h);
guidata(hObject,h);

function h = skewFig(h)
hval0            = [h.dat.stat.skew];
hval             = hval0 / nanmean(hval0);
hval             = max(0, hval - (1 - 2*nanstd(hval))) + 1;
hval             = log(hval) / nanmean(log(hval));
hval             = .6 * (1 - min(2*nanstd(hval)+1, hval)/(nanstd(hval)*2 + 1)) + .15;
h.dat.cl.rands   = hval;
h.dat.cl.cmap    = [hval0(:) hval(:)]; 
I = redraw_figure(h);
h.dat.cl.skewFig = I;
set_maskCcolor(h, 3);

% --- ellipse
function pushbutton112_Callback(hObject, eventdata, h)
h=ellipseFig(h);
guidata(hObject,h);

function h=ellipseFig(h)
hval0            = min(5, [h.dat.stat.aspect_ratio]);
hval             = hval0;
hval             = hval - min(hval);
hval             = hval / nanmean(hval);
hval             = .6 * (1 - min(2*nanstd(hval)+1, hval)/(nanstd(hval)*2 + 1)) + .15;
h.dat.cl.rands   = hval;
h.dat.cl.cmap    = [hval0(:) hval(:)]; 
I = redraw_figure(h);
h.dat.cl.ellipseFig = I;
set_maskCcolor(h, 8);

% --- meanimg
function pushbutton102_Callback(hObject, eventdata, h)
h=meanimgFig(h);
guidata(hObject,h);

function h=meanimgFig(h)
hval0            = [h.dat.stat.mimgProjAbs];
hval             = hval0;
hval             = hval - min(hval);
hval             = hval / nanmean(hval);
hval             = .6 * (1 - min(2*nanstd(hval)+1, hval)/(nanstd(hval)*2 + 1)) + .15;
h.dat.cl.rands   = hval;
h.dat.cl.cmap    = [hval0(:) hval(:)]; 
I = redraw_figure(h);
h.dat.cl.meanimgFig = I;
set_maskCcolor(h, 4);


% --- cmpct
function pushbutton99_Callback(hObject, eventdata, h)
h=cmpctFig(h);
guidata(hObject,h);

function h=cmpctFig(h)
hval0            = min(3, [h.dat.stat.cmpct]);
hval             = hval0;
hval             = hval - min(hval);
hval             = hval / nanmean(hval);
hval             = .6 * (1 - min(2*nanstd(hval)+1, hval)/(nanstd(hval)*2 + 1)) + .15;
h.dat.cl.rands   = hval;
h.dat.cl.cmap    = [hval0(:) hval(:)]; 
I = redraw_figure(h);
h.dat.cl.cmpctFig = I;
set_maskCcolor(h, 5);

% --- footprint
function pushbutton100_Callback(hObject, eventdata, h)
h=footprintFig(h);
guidata(hObject,h);

function h=footprintFig(h)
hval0            = [h.dat.stat.footprint];
hval             = hval0;
hval             = hval - min(hval);
hval             = hval / nanmean(hval);
hval             = .6 * (1 - min(2*nanstd(hval)+1, hval)/(nanstd(hval)*2 + 1)) + .15;
h.dat.cl.rands   = hval;
h.dat.cl.cmap    = [hval0(:) hval(:)]; 
I = redraw_figure(h);
h.dat.cl.footprintFig = I;
set_maskCcolor(h, 6);


% --- red channel
function pushbutton104_Callback(hObject, eventdata, h)
h=redFig(h);
guidata(hObject,h);

function h=redFig(h)
hval0            = [h.dat.stat.redprob];
hval             = hval0;
hval             = hval / nanmean(hval);
hval             = max(0, hval - (1 - 2*nanstd(hval))) + 1;
hval             = log(hval) / nanmean(log(hval));
hval             = .6 * (1 - min(2*nanstd(hval)+1, hval)/(nanstd(hval)*2 + 1)) + .15;
h.dat.cl.rands   = hval;
h.dat.cl.cmap    = [hval0(:) hval(:)]; 
I = redraw_figure(h);
h.dat.cl.redFig  = I;
set_maskCcolor(h, 7);



% ------------------ CLASSIFIER --------------------------%

% choose classifier? button
function pushbutton93_Callback(hObject, eventdata, h)
rootS2p = which('run_pipeline');
if isempty(rootS2p)
    error('could not identify Suite2p location! where is new_main.m?')
end
rootS2p = fileparts(rootS2p);
rootS2p = fullfile(rootS2p, 'configFiles');

[filename1,filepath1]   = uigetfile(fullfile(rootS2p, '*.mat'), 'Select classifier file');
if filename1
    h.dat.cl.fpath          = fullfile(filepath1, filename1);
    h                       = classROI(h);
    
    redraw_figure(h);
    
%     hload = load(h.dat.cl.fpath);
%     h.st = hload.st;
%     h.save_cls = h.dat.cl.fpath;
%     h.prior = prior;

    guidata(hObject,h);
end

% --- threshold setting for classifier
function edit50_Callback(hObject, eventdata, h)
h.dat.cl.threshold = str2double(get(h.edit50,'String'));
for j = 1:length(h.dat.stat)
    h.dat.stat(j).iscell = h.dat.stat(j).cellProb > h.dat.cl.threshold;
end
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

% --- new classifier button
function pushbutton101_Callback(hObject, eventdata, h)
rootS2p = which('run_pipeline');
if isempty(rootS2p)
    warndlg('Suite2p location not in PATH! where is run_pipeline.m?')
    [filename1,rootS2p]=uigetfile(root, 'Select run_pipeline.m');
    
end
rootS2p = fileparts(rootS2p);
rootS2p = fullfile(rootS2p, 'configFiles');

def_name = fullfile(rootS2p, 'cl_new.mat');

[filename1,filepath1]   = uigetfile(fullfile(rootS2p, 'priors', '*.mat'), ...
        'Select priors file, will be used to seed the new classifier');
if filename1
    prior_file = fullfile(filepath1, filename1);
else
    error('you must choose a priors file')
end

[FileName,PathName] = uiputfile('*.mat', 'Create new classifier', def_name); 

if FileName
    load(prior_file);
    st                      = [];
    save(fullfile(PathName, FileName), 'st', 'prior', 'statLabels')
    
    h.dat.cl.fpath          = fullfile(PathName, FileName);
    h                       = classROI(h);
    
    redraw_figure(h);

    hload = load(h.dat.cl.fpath);
    h.st = hload.st;
    h.save_cls = h.dat.cl.fpath;
    
    guidata(hObject,h);
else
    warndlg('you did not create a new classifier file');
    error('you did not create a new classifier file')
end

%----------------------------- HELP SECTION ----------------------%

function pushbutton64_Callback(hObject, eventdata, h)

msg{1} = ['First, you need to load an F_*_*_*.mat produced by Suite2p. ']; 
msg{3} = ['Now start exploring the data with the various visualizations. Press the help buttons next to each category to learn the options.'];
msg{5} = ['You can label an ROI as cell/not cell, which determines the image axes it is shown in.'];
msg{7} = ['Right-click on an ROI to flip its label.'];
msg{9} = ['The more you do this, the better the automated classification becomes.'];
msg{11} = ['hint: you can keep any help box open while interacting with the GUI. '];

msgbox(msg, 'Instructions!');

function pushbutton108_Callback(hObject, eventdata, handles)
msg{1} = ['Timecourse of selected ROI and neuropil fluorescence across experiment.'];
msg{3} = ['The selected cell has 0 saturation (gray) in ROI image.'];
msg{5} = ['Traces have been smoothed for display purposes.'];
msg{7} = ['Statistics shown are same variables as used in "mask color" section.'];

msgbox(msg, 'Fluorescence instructions');

function pushbutton106_Callback(hObject, eventdata, handles)
msg{1} = ['This applies only if you select "ROIs" under "background".'];
msg{3} = ['RANDOM: color chosen randomly - RED cells labelled in red (if secondary channel)'];
msg{5} = ['*** middle mouse-click on an ROI to switch RED label ON/OFF ***'];
msg{7} = ['for all other selections, color/hue for each ROI varies purple to yellow'];
msg{8} = ['(scale bar for colors shown below buttons):'];
msg{10} = ['CLASSIFIER: probability assigned by classifier'];
msg{11} = ['SKEW: skewness of activity, after neuropil correction and some smoothing'];
msg{12} = ['MEANIMG: weighting of activity mask onto mean image'];
msg{13} = ['CMPCT: compactness of ROI pixels. Smallest is 1, for disks.'];
msg{14} = ['FOOT: "footprint" of ROI; ~ number of correlated neighboring pixels'];
msg{15} = ['RED: probability of being a red-tagged cell, assigned by algorithm'];
msg{17} = ['hint: the letters in paranthesis are keyboard shortcuts.'];

msgbox(msg, 'Mask color instructions');

function pushbutton107_Callback(hObject, eventdata, handles)
msg{1} = ['This applies only if you select "ROIs" under "background".'];
msg{3} = ['Selection determines the brightness for all pixels inside ROIs.'];
msg{5} = ['UNIT NORM: values are normalized to unit norm per ROI'];
msg{7} = ['VARIANCE: fraction of variance explained in each pixel by parent ROI'];

msgbox(msg, 'Mask brightness instructions');

function pushbutton105_Callback(hObject, eventdata, handles)
msg{1} = ['ROI: shows the masks identified by Suite2p, colored according to the property selected under "mask color"'];
msg{3} = ['CORR: shows the correlation map between a pixel and its nearby pixels'];
msg{5} = ['MEAN: average registered image'];
msg{7} = ['RED: average registered image of "red"/secondary color channel'];
msg{9} = ['RED - GREEN: subtracts off the contamination from "green"/primary channel'];
msg{11} = ['PROC: switch to toggle image contrast normalization'];
msg{13} = ['hint: the letters in paranthesis are keyboard shortcuts'];

msgbox(msg, 'Background instructions');


function pushbutton110_Callback(hObject, eventdata, handles)
msg{1} = ['The classifier assigns the initial "iscell" labels to the ROIs.'];
msg{3} = ['Initially, it might not work very well, but will improve as you make choices in the GUI (and save the "proc" files).'];
msg{5} = ['To begin, press "new classifier?" and choose one of the provided priors.'];
msg{7} = ['The last used classifier will be automatically selected, every time you load a new plane.'];
msg{9} = ['As you process datasets, a database of your cells is built, and the classifier is re-trained.'];
msg{11} = ['The prior counts for about 300 cells, so your choices will start making a difference after a few hundred manually validated cells.'];
msg{13} = ['The system allows you to build personalized classifiers for different kinds of data (i.e. somas, dendrites, boutons, different brain areas, calcium indicators, or zoom levels).'];

msgbox(msg, 'Classifier instructions');

function pushbutton111_Callback(hObject, eventdata, handles)
msg{1} = ['Zoom in on portion of the image. Quadrants have 10% overlap.'];

msg{3} = ['Buttons become dark grey after visiting a quadrant.'];

msgbox(msg, 'ZOOM panel instructions');




% --- Executes on selection change in popupmenu12.
function popupmenu12_Callback(hObject, eventdata, h)
% hObject    handle to popupmenu12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(hObject,'String'));
popstr = contents{get(hObject,'Value')};

switch popstr
    case 'number of pixels'
        stn = 1;
    case 'compactness'
        stn = 2;
    case 'aspect ratio'
        stn = 3;
    case 'skewness'
        stn = 4;
    case 'standard dev'
        stn = 5;
    case 'footprint'
        stn = 6;
    case 'meanimg proj'
        stn = 7;
end
set(h.edit52,'String',num2str(h.statmins(stn)));
set(h.edit54,'String',num2str(h.statmaxs(stn)));

h.statstr = h.statstrs{stn};
h.statnum = stn;

guidata(hObject,h);

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu12 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu12


% --- Executes during object creation, after setting all properties.
function popupmenu12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit52_Callback(hObject, eventdata, h)
% hObject    handle to edit52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h.statmins(h.statnum) = str2double(get(hObject,'String'));

goodcells = set_thres(h.dat.stat, h.statstrs, h.statmins, h.statmaxs);
[h.dat.stat(~goodcells).iscell] = deal(0);
[h.dat.stat(goodcells).iscell]  = deal(1);

redraw_figure(h);

guidata(hObject,h);
% Hints: get(hObject,'String') returns contents of edit52 as text
%        str2double(get(hObject,'String')) returns contents of edit52 as a double


% --- Executes during object creation, after setting all properties.
function edit52_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit54_Callback(hObject, eventdata, h)
% hObject    handle to edit54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h.statmaxs(h.statnum) = str2double(get(hObject,'String'));

goodcells = set_thres(h.dat.stat, h.statstrs, h.statmins, h.statmaxs);
[h.dat.stat(~goodcells).iscell] = deal(0);
[h.dat.stat(goodcells).iscell]  = deal(1);

redraw_figure(h);

guidata(hObject,h);

% Hints: get(hObject,'String') returns contents of edit54 as text
%        str2double(get(hObject,'String')) returns contents of edit54 as a double


% --- Executes during object creation, after setting all properties.
function edit54_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function goodcells = set_thres(stat, statstrs, statmins, statmaxs)
goodcells = true(size(stat));
for j = 1:length(statstrs) 
    svals = [stat.(statstrs{j})];
    goodcells = goodcells & (svals > statmins(j) & svals < statmaxs(j));    
end
    
    
    




