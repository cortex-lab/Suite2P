function varargout = registers2p(varargin)
% REGISTERS2P MATLAB code for registers2p.fig
%      REGISTERS2P, by itself, creates a new REGISTERS2P or raises the existing
%      singleton*.
%
%      H = REGISTERS2P returns the handle to a new REGISTERS2P or the handle to
%      the existing singleton*.
%
%      REGISTERS2P('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REGISTERS2P.M with the given input arguments.
%
%      REGISTERS2P('Property','Value',...) creates a new REGISTERS2P or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before registers2p_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to registers2p_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help registers2p

% Last Modified by GUIDE v2.5 26-May-2017 15:46:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @registers2p_OpeningFcn, ...
                   'gui_OutputFcn',  @registers2p_OutputFcn, ...
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


% --- Executes just before registers2p is made visible.
function registers2p_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to registers2p (see VARARGIN)

% Choose default command line output for registers2p
handles.output = hObject;
num_ims = 2;
handles.image_display = [];
handles.roi_display = [];
handles.cursors(1).scatter = [];
handles.cursors(2).scatter = [];
handles.targets(1).scatter = [];
handles.targets(2).scatter = [];
handles.dil_targets_handles = [];
handles.original_targets = [];
handles.plotted_targets = [];
handles.overlay = false;
handles.im_flag = true;
handles.colors = [0 0 0];
handles.target_sizes = [15 70];
handles.bg_color = [1 1 1 ; 0 1 0];
handles.mask_fields = fields(handles);
handles.mask_fields = handles.mask_fields(cellfun(@(y) ~isempty(strfind(y,'mask_text')),handles.mask_fields));
handles.mask_field_nums = cellfun(@(x) str2num(x(strfind(x,'text')+4:end)),handles.mask_fields);
for i = 1:num_ims
    handles.files(i).filepath = [];
    % Variables
    handles.files(i).num_masks = 5;
    handles.files(i).image{1} = zeros(512,512,handles.files(i).num_masks);
    handles.files(i).image{2} = handles.files(i).image{1};
    handles.files(i).dims = size(handles.files(i).image{1});
    handles.files(i).dims = handles.files(i).dims(1:2);
    handles.files(i).rois_ibar = [];
    handles.files(i).centroids = [];
    handles.files(i).iscell_idcs = [];
    % Calculated
    handles.files(i).roi_status = [];
    handles.files(i).roi_idcs = [];
    handles.files(i).num_rois = [];
    handles.files(i).roi_areas = [];
    handles.files(i).iscell = [];
    handles.files(i).loaded = false;
    handles.files(i).targets_overlap_yx = [nan nan];
    % Defaults
    handles.files(i).default_image = handles.files(i).image;
    handles.files(i).default_dims = [];
    handles.files(i).default_rois_ibar = [];
    handles.files(i).default_centroids = [];
    handles.files(i).default_roi_idcs = [];
    handles.files(i).default_roi_areas = [];
    handles.files(i).num_rois = [];
end
handles.overlap_thold = 0.6;
handles.targets_overlap_thold = 5;
handles = set_gui_defaults(hObject,handles);

handles.transform_filepath = [];
handles.tform = [];
handles.tform_direction = [0 1];
handles.transform_loaded = false;

handles.targets_filepath = [];
handles.targets_image = [];
handles.dil_targets_yx = [];
handles.tformed_targets_yx = [nan nan];
handles.tformed_targets_default = [nan nan];
handles.targets_loaded = false;
handles.targets_s2p_rois = [];
handles.selected_targets = false;
handles.targets_overlapped = false;

handles.displayed_image = 1;
handles.displayed_mask = 1;
handles.proc_flag = 1;
handles.image_display = imagesc(handles.image_axis,handles.files(1).image{1}(:,:,1),'ButtonDownFcn',@lasso_targets);
colormap(handles.image_axis,'gray')
axis(handles.image_axis,'square','off')
hold(handles.image_axis,'on');
handles.image_axis.Units = 'normalized';
handles.targets(1).scatter = scatter(handles.image_axis,[],[],30,'r.');
handles.targets(2).scatter = scatter(handles.image_axis,[],[],30,'r.');
handles.original_targets = scatter(handles.image_axis,[],[],20,[0.6 0.6 0.6],'o');
handles.plotted_targets = scatter(handles.image_axis,[],[],handles.target_sizes(1),[1 0 1],'o','filled');
handles.cursors(1).scatter = scatter(handles.image_axis,nan,nan,50,'o','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[0 0.5 1]);

handles.roi_display = image(handles.roi_axis,handles.files(1).image{1}(:,:,1),'ButtonDownFcn',@select_rois);
colormap(handles.roi_axis,'gray')
axis(handles.roi_axis,'square','off')
hold(handles.roi_axis,'on');
handles.image_axis.Units = 'normalized';
handles.cursors(2).scatter = scatter(handles.roi_axis,nan,nan,50,'o','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[0 0.5 1]);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes registers2p wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = registers2p_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function fixed_textbox_Callback(hObject, eventdata, handles)
% hObject    handle to fixed_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fixed_textbox as text
%        str2double(get(hObject,'String')) returns contents of fixed_textbox as a double


% --- Executes during object creation, after setting all properties.
function fixed_textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fixed_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function transform_textbox_Callback(hObject, eventdata, handles)
% hObject    handle to transform_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of transform_textbox as text
%        str2double(get(hObject,'String')) returns contents of transform_textbox as a double


% --- Executes during object creation, after setting all properties.
function transform_textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to transform_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function moving_textbox_Callback(hObject, eventdata, handles)
% hObject    handle to moving_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of moving_textbox as text
%        str2double(get(hObject,'String')) returns contents of moving_textbox as a double


% --- Executes during object creation, after setting all properties.
function moving_textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to moving_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in load_fixed_button.
function load_fixed_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_fixed_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
idx = 1;
handles = load_s2p_struct(hObject,handles,idx);
guidata(hObject,handles)


% --- Executes on button press in load_moving_button.
function load_moving_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_moving_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
idx = 2;
handles = load_s2p_struct(hObject,handles,idx);
guidata(hObject,handles)



% --- Executes on button press in load_transform_button.
function load_transform_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_transform_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.files(1).filepath) && ~isempty(handles.files(2).filepath)
    name_idcs = find(cellfun(@(x) ~isempty(x),{handles.files(:).filepath}));
    dd = fileparts(handles.files(name_idcs(1)).filepath);
    try
        [f,d] = uigetfile([dd filesep '*.mat']);
    catch
        [f,d] = uigetfile('*.mat');
    end
    if f~=0
        handles.transform_filepath = [d filesep f];
        handles.transform_textbox.String = handles.transform_filepath;
        load(handles.transform_filepath)
        handles.tform = t.tform;
        handles.tform_direction = t.tform_direction;
        handles.transform_loaded = true;
        
        order = handles.tform_direction+1;
        handles = set_gui_defaults(hObject,handles);
        handles = reload_roi_defaults(handles);
        i = 0;
        h = waitbar(i,'Fitting transform...');
        steps = (handles.files(order(2)).num_masks * 2) + 3;
        i = i+1;
        waitbar(i / steps,h,'Registering images...')
        for p = 1:2
            for m = 1:handles.files(order(2)).num_masks
                handles.files(order(2)).image{p}(:,:,m) = imwarp(handles.files(order(2)).image{p}(:,:,m),handles.tform,'OutputView',imref2d(size(handles.files(order(2)).image{1}(:,:,1))));
                i = i+1;
                waitbar(i / steps);
            end
        end
        i = i+1;
        waitbar(i / steps,h,'Registering ROIs...')
        [handles] = transform_rois(handles,order(2));
        if handles.targets_loaded && handles.targets_overlapped
            handles = reset_targets_overlap(handles);
        end
        waitbar(1);
        close(h);
        update_displays(hObject,handles);
        guidata(hObject,handles)
    end
end


% --- Executes on button press in load_targets_button.
function load_targets_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_targets_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[f,d] = uigetfile('*.*');
if f~=0
    handles.targets_filepath = [d f];
    handles.targets_textbox.String = handles.targets_filepath;
    [~,~,ext] = fileparts(handles.targets_filepath);
    switch ext
        case {'.tif' '.tiff'}
            handles.targets_image = imread(handles.targets_filepath);
            handles.tformed_targets_yx = [];
            [handles.tformed_targets_yx(:,1),handles.tformed_targets_yx(:,2)] = find(handles.targets_image);
            handles.tformed_targets_default = handles.tformed_targets_yx;
            handles.targets_loaded = true;
            handles.targets_s2p_rois = zeros(size(handles.tformed_targets_yx,1),3);
            handles.colors = repmat([1 0 0],[size(handles.tformed_targets_yx(:,1)),1]);
            handles.selected_targets = false * handles.tformed_targets_yx(:,1);
        case {'.mat'}
            load(handles.targets_filepath)
            if size(params.targets_yxzc,2) == 4
                [~,order] = sort(params.targets_yxzc(:,4));
            else
                order = 1:size(params.targets_yxzc,1);
            end
            handles.tformed_targets_yx = params.targets_yxzc(order,1:2);
            handles.tformed_targets_default = handles.tformed_targets_yx;
            handles.targets_loaded = true;
            handles.targets_s2p_rois = zeros(size(handles.tformed_targets_yx,1),3);
            handles.colors = repmat([1 0 0],[size(handles.tformed_targets_yx(:,1)),1]);
            handles.selected_targets = false * handles.tformed_targets_yx(:,1);
    end
    if handles.targets_loaded
        if handles.targets_overlapped
            handles = reset_targets_overlap(handles);
        end
        handles = plot_targets(hObject,handles);
    end
end
guidata(hObject,handles)


% --- Executes on button press in save_tform_button.
function save_tform_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_tform_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.files(1).filepath) && ~isempty(handles.files(2).filepath) && ~isempty(handles.tform)
    name_idcs = find(cellfun(@(x) ~isempty(x),{handles.files(:).filepath}));
    dd = fileparts(handles.files(name_idcs(1)).filepath);
    name_idcs = find(cellfun(@(x) ~isempty(x),{handles.files(:).filepath}));
    name_chunks = cell(1,numel(name_idcs));
    ff = name_chunks;
    for i = 1:numel(name_idcs)
        [dd,ff{i}] = fileparts(handles.files(name_idcs(i)).filepath);
        name_chunks{i} = ff{i}(3:end);
        name_chunks{i} = strrep(name_chunks{i},'proc','');
    end
    default_name = [name_chunks{:}];
    default_name(end) = [];
    default_name = [default_name '_tform'];
    [f,d] = uiputfile([dd filesep '*.mat'],'Save...',[dd filesep default_name '.mat']);
    if d ~= 0
        t.tform = handles.tform;
        t.tform_direction = handles.tform_direction;
        save([d filesep f],'t');
    end
end



% --- Executes on button press in register_button.
function register_button_Callback(hObject, eventdata, handles)
% hObject    handle to register_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.files(1).loaded && handles.files(2).loaded 
    order = handles.tform_direction+1;
    moving_im  = double(handles.files(order(2)).default_image{handles.proc_flag}(:,:,handles.displayed_mask));
    static_im = double(handles.files(order(1)).default_image{handles.proc_flag}(:,:,handles.displayed_mask));
    moving_im = moving_im - min(moving_im(:));
    static_im = static_im - min(static_im(:));
    moving_im  = uint8(moving_im/max(max(moving_im))*255);
    static_im = uint8(static_im/max(max(static_im))*255);
    
    [movingPoints, fixedPoints] = cpselect(moving_im, static_im, 'wait',true);
    if ~isempty(movingPoints) && ~isempty(fixedPoints)
        i = 0;
        h = waitbar(i,'Fitting transform...');
        handles.tform = fitgeotrans(movingPoints, fixedPoints, 'affine');
        handles = set_gui_defaults(hObject,handles);
        handles = reload_roi_defaults(handles);
        steps = (handles.files(order(2)).num_masks * 2) + 3;
        i = i+1;
        waitbar(i / steps,h,'Registering images...')
        for p = 1:2
            for m = 1:handles.files(order(2)).num_masks
                handles.files(order(2)).image{p}(:,:,m) = imwarp(handles.files(order(2)).image{p}(:,:,m),handles.tform,'OutputView',imref2d(size(handles.files(order(2)).image{1}(:,:,1))));
                i = i+1;
                waitbar(i / steps);
            end
        end
        i = i+1;
        waitbar(i / steps,h,'Registering ROIs...')
        handles = transform_rois(handles,order(2));
        if handles.targets_loaded && handles.targets_overlapped
            handles = reset_targets_overlap(handles);
        end
        update_displays(hObject,handles);
        waitbar(1);
        close(h);
    end
end
guidata(hObject,handles);


% --- Executes on button press in reset_button.
function reset_button_Callback(hObject, eventdata, handles)
% hObject    handle to reset_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = reload_roi_defaults(handles);   

handles.transform_filepath = [];
handles.tform = [];
handles.tform_direction = [0 1];
handles.transform_loaded = false;
handles.transform_textbox.String = '...';

handles.tformed_targets_yx = handles.tformed_targets_default;
handles = reset_rois_overlap(handles);
if handles.targets_loaded && handles.targets_overlapped
    handles = reset_targets_overlap(handles);
end
update_displays(hObject,handles);
guidata(hObject,handles);


function targets_textbox_Callback(hObject, eventdata, handles)
% hObject    handle to targets_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of targets_textbox as text
%        str2double(get(hObject,'String')) returns contents of targets_textbox as a double


% --- Executes during object creation, after setting all properties.
function targets_textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to targets_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function image_axis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate image_axis


% --- Executes on mouse press over axes background.
function roi_axis_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to roi_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on mouse press over axes background.
function image_axis_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to image_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.files(1).filepath) || ~isempty(handles.files(2).filepath)
    name_idcs = find(cellfun(@(x) ~isempty(x),{handles.files(:).filepath}));
    name_chunks = cell(1,numel(name_idcs));
    ff = name_chunks;
    for i = 1:numel(name_idcs)
        [dd,ff{i}] = fileparts(handles.files(name_idcs(i)).filepath);
        name_chunks{i} = ff{i}(3:end);
        name_chunks{i} = strrep(name_chunks{i},'proc','');
    end
    default_name = [name_chunks{:}];
    default_name(end) = [];
    default_name = [default_name '_reg'];
    [f,d] = uiputfile([dd filesep '*.mat'],'Save...',[dd filesep default_name '.mat']);
    if d ~= 0
        fnames = fieldnames(handles);
        s_idx = find(cellfun(@(x) strcmp(x,'files'),fnames));
        regi = [];
        for fn = s_idx:numel(fnames)
            regi.dat.(fnames{fn}) = handles.(fnames{fn});
        end
        regi.dat.overlapping_rois(1,:) = [];
        regi.dat.p_overlaps(1) = [];
        regi.rois.idcs = regi.dat.overlapping_rois;
        regi.targets.idcs = regi.dat.targets_s2p_rois;
        regi.rois.iscell_idcs = 0 * regi.dat.overlapping_rois;
        regi.targets.iscell_idcs = nan * regi.dat.targets_s2p_rois;
        for i = 1:numel(regi.dat.files)
            for r = 1:size(regi.dat.overlapping_rois(:,i))
                regi.rois.iscell_idcs(r,i) = find(handles.files(i).iscell == regi.dat.overlapping_rois(r,i));
            end
            if handles.targets_loaded
                for t = 1:size(regi.dat.targets_s2p_rois(:,i))
                    if ~isnan(regi.dat.targets_s2p_rois(t,i))
                        regi.targets.iscell_idcs(t,i) = find(handles.files(i).iscell == regi.dat.targets_s2p_rois(t,i));
                    end
                end
            end
        end
        if handles.targets_loaded
            regi.targets.iscell_idcs(:,3) = regi.dat.targets_s2p_rois(:,3);
        end
        save([d filesep f],'regi');
    end
else
    uiwait(msgbox('No data loaded - please load at least 1 session before trying to save', 'Error','error'));
end


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

if isstrprop(eventdata.Key,'digit')
    num = str2num(eventdata.Key);
    if max([1:handles.files(handles.displayed_image).num_masks] == num) == 1
        handles.displayed_mask = num;
         handles.(handles.mask_fields{handles.mask_field_nums == num}).BackgroundColor = handles.bg_color(2,:);
         off = handles.mask_fields(handles.mask_field_nums ~= num);
         for i = 1:numel(off)
             handles.(off{i}).BackgroundColor = handles.bg_color(1,:);
         end
        update_displays(hObject,handles);
    end
else
    switch eventdata.Key
        %     case '1'
        %         handles.displayed_mask = 2;
        %         update_displays(hObject,handles);
        %     case '2'
        %         handles.displayed_mask = 5;
        %         update_displays(hObject,handles);
        case 'p'
            handles.proc_flag = ~(handles.proc_flag-1)+1;
            handles.proc_text.BackgroundColor = handles.bg_color(handles.proc_flag,:);
            update_displays(hObject,handles);
        case 'q'
            handles.im_flag = true;
            handles.displayed_image = 1;
            handles.ROI_text1.BackgroundColor = handles.bg_color(1,:);
            handles.ROI_text2.BackgroundColor = handles.bg_color(1,:);
            handles.day_text1.BackgroundColor = handles.bg_color(2,:);
            handles.day_text2.BackgroundColor = handles.bg_color(1,:);
            update_displays(hObject,handles);
        case 'w'
            handles.im_flag = true;
            handles.displayed_image = 2;
            handles.ROI_text1.BackgroundColor = handles.bg_color(1,:);
            handles.ROI_text2.BackgroundColor = handles.bg_color(1,:);
            handles.day_text1.BackgroundColor = handles.bg_color(1,:);
            handles.day_text2.BackgroundColor = handles.bg_color(2,:);
            update_displays(hObject,handles);
        case 'a'
            handles.im_flag = false;
            handles.displayed_image = 1;
            handles.ROI_text1.BackgroundColor = handles.bg_color(2,:);
            handles.ROI_text2.BackgroundColor = handles.bg_color(1,:);
            handles.day_text1.BackgroundColor = handles.bg_color(1,:);
            handles.day_text2.BackgroundColor = handles.bg_color(1,:);
            update_displays(hObject,handles);
        case 's'
            handles.im_flag = false;
            handles.displayed_image = 2;
            handles.ROI_text1.BackgroundColor = handles.bg_color(1,:);
            handles.ROI_text2.BackgroundColor = handles.bg_color(2,:);
            handles.day_text1.BackgroundColor = handles.bg_color(1,:);
            handles.day_text2.BackgroundColor = handles.bg_color(1,:);
            update_displays(hObject,handles);
        case 'o'
            handles.overlay = ~handles.overlay;
            handles.cursor_overlay_text.BackgroundColor = handles.bg_color(handles.overlay+1,:);
            guidata(hObject,handles)
            if handles.overlay
                order = [1 0];
                C = get (handles.image_axis, 'CurrentPoint');
                C = C(1,1:2);
                if (C(1) > handles.files(1).dims(1) | C(2) > handles.files(1).dims(2)) || min(C)<0
                    handles.cursors(2).scatter.XData = nan;
                    handles.cursors(2).scatter.YData = nan;
                    C = get (handles.roi_axis, 'CurrentPoint');
                    C = C(1,1:2);
                    idcs = (~order)+1;
                else
                    idcs = order+1;
                end
                if (C(1) <= handles.files(idcs(1)).dims(1) | C(2) <= handles.files(idcs(1)).dims(2)) || min(C)>0
                    handles.cursors(idcs(1)).scatter.XData = C(1,1);
                    handles.cursors(idcs(1)).scatter.YData = C(1,2);
                else
                    handles.cursors(idcs(1)).scatter.XData = nan;
                    handles.cursors(idcs(1)).scatter.YData = nan;
                end
            else
                for i = 1:numel(handles.cursors)
                    handles.cursors(i).scatter.XData = nan;
                    handles.cursors(i).scatter.YData = nan;
                end
            end
    end
    
    if ~isempty(handles.tformed_targets_yx) && max(handles.selected_targets) == 1
        switch eventdata.Key
            case 'downarrow'
                if max(handles.tformed_targets_yx(handles.selected_targets,1)+1) < handles.files(handles.displayed_image).dims(1) %513
                    handles.tformed_targets_yx(handles.selected_targets,1) = handles.tformed_targets_yx(handles.selected_targets,1)+1;
                end
                handles.plotted_targets.XData(handles.selected_targets) = handles.tformed_targets_yx(handles.selected_targets,2);
                handles.plotted_targets.YData(handles.selected_targets) = handles.tformed_targets_yx(handles.selected_targets,1);
            case 'uparrow'
                if min(handles.tformed_targets_yx(handles.selected_targets,1)-1) > 0
                    handles.tformed_targets_yx(handles.selected_targets,1) = handles.tformed_targets_yx(handles.selected_targets,1)-1;
                end
                handles.plotted_targets.XData(handles.selected_targets) = handles.tformed_targets_yx(handles.selected_targets,2);
                handles.plotted_targets.YData(handles.selected_targets) = handles.tformed_targets_yx(handles.selected_targets,1);
            case 'leftarrow'
                if min(handles.tformed_targets_yx(handles.selected_targets,2)-1) > 0
                    handles.tformed_targets_yx(handles.selected_targets,2) = handles.tformed_targets_yx(handles.selected_targets,2)-1;
                end
                handles.plotted_targets.XData(handles.selected_targets) = handles.tformed_targets_yx(handles.selected_targets,2);
                handles.plotted_targets.YData(handles.selected_targets) = handles.tformed_targets_yx(handles.selected_targets,1);
            case 'rightarrow'
                if max(handles.tformed_targets_yx(handles.selected_targets,2)+1) < handles.files(handles.displayed_image).dims(2) %513
                    handles.tformed_targets_yx(handles.selected_targets,2) = handles.tformed_targets_yx(handles.selected_targets,2)+1;
                end
                handles.plotted_targets.XData(handles.selected_targets) = handles.tformed_targets_yx(handles.selected_targets,2);
                handles.plotted_targets.YData(handles.selected_targets) = handles.tformed_targets_yx(handles.selected_targets,1);
            case {'escape' 'return'}
                handles.selected_targets = false(numel(handles.selected_targets),1);
                handles.plotted_targets.CData(:,2) = handles.selected_targets;
            case 'd'
                handles.tformed_targets_yx(handles.selected_targets,:) = handles.tformed_targets_default(handles.selected_targets,:);
                handles.plotted_targets.XData(handles.selected_targets) = handles.tformed_targets_yx(handles.selected_targets,2);
                handles.plotted_targets.YData(handles.selected_targets) = handles.tformed_targets_yx(handles.selected_targets,1);
        end
        if handles.targets_overlapped
            handles = detect_target_overlap(handles);
        end
    end
end
guidata(hObject,handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CustomFunctions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = load_s2p_struct(hObject,handles,idx)
[f,d] = uigetfile('*.*');
if f~=0
    handles = set_gui_defaults(hObject,handles);
    handles.files(idx).filepath = [d f];
    handles.fixed_textbox.String = handles.files(1).filepath;
    handles.moving_textbox.String = handles.files(2).filepath;
    [~,~,ext] = fileparts(handles.files(idx).filepath);
    switch ext
        case '.mat'
            try
                load(handles.files(idx).filepath);

                handles.files(idx).dims = size(dat.ops.mimg1);
                handles.files(idx).dims = handles.files(idx).dims(1:2);
                %handles.files(idx).num_masks = size(dat.mimg,3);
                handles.files(idx).num_masks = 5;
                handles.files(idx).image{1} = zeros(handles.files(idx).dims(1),handles.files(idx).dims(2),handles.files(idx).num_masks);
                handles.files(idx).image{2} = zeros(handles.files(idx).dims(1),handles.files(idx).dims(2),handles.files(idx).num_masks);
                
                [xx,yy] = meshgrid(dat.ops.xrange,dat.ops.yrange);
                im_idcs = sub2ind(size(handles.files(idx).image{1}(:,:,1)),yy,xx);
                pixels_per_im = handles.files(idx).dims(1) * handles.files(idx).dims(2);
                
%                 for i = 1:handles.files(idx).num_masks
%                     handles.files(idx).image{1}(im_idcs + ((i-1)*pixels_per_im)) = dat.mimg(:,:,i);
%                     handles.files(idx).image{2}(im_idcs + ((i-1)*pixels_per_im)) = dat.mimg_proc(:,:,i);
%                 end
                
                % ROI_text1 im
                i = 1;
                handles.files(idx).image{1}(im_idcs + ((i-1)*pixels_per_im)) = dat.mimg(:,:,2);
                handles.files(idx).image{2}(im_idcs + ((i-1)*pixels_per_im)) = dat.mimg_proc(:,:,2);
                % corr im
                i = 2;
                handles.files(idx).image{1}(im_idcs + ((i-1)*pixels_per_im)) = dat.mimg(:,:,5);
                handles.files(idx).image{2}(im_idcs + ((i-1)*pixels_per_im)) = dat.mimg_proc(:,:,5);
                % sd im
                i = 3;
                handles.files(idx).image{1}(im_idcs + ((i-1)*pixels_per_im)) = dat.ops.sdmov;
                handles.files(idx).image{2}(im_idcs + ((i-1)*pixels_per_im)) = dat.ops.sdmov;
                % corr im
                i = 4;
                handles.files(idx).image{1}(im_idcs + ((i-1)*pixels_per_im)) = dat.res.iclust;
                handles.files(idx).image{2}(im_idcs + ((i-1)*pixels_per_im)) = dat.res.iclust;
                % skewness im
                i = 5;
                handles.files(idx).image{1}(im_idcs + ((i-1)*pixels_per_im)) = dat.res.lambda;
                handles.files(idx).image{2}(im_idcs + ((i-1)*pixels_per_im)) = dat.res.lambda0;
                
                handles.files(idx).default_image = handles.files(idx).image;
                
                handles = set_gui_defaults(hObject,handles);
                handles = reload_roi_defaults(handles);
                if handles.targets_loaded && handles.targets_overlapped
                    handles = reset_targets_overlap(handles);
                end
                
                [handles.files(idx).rois_ibar,handles.files(idx).centroids,handles.files(1).roi_areas,handles.files(idx).iscell] = load_S2P_rois(dat,handles.files(idx).dims);
                handles.files(idx).roi_idcs = unique(handles.files(idx).rois_ibar(:,4),'stable');
                handles.files(idx).num_rois = numel(unique(handles.files(idx).rois_ibar(:,4)));
                handles.files(idx).roi_status = [ones(size(handles.files(idx).rois_ibar,1),1) zeros(size(handles.files(idx).rois_ibar,1),1)];
                
                handles.files(idx).default_rois_ibar = handles.files(idx).rois_ibar;
                handles.files(idx).default_centroids = handles.files(idx).centroids;
                handles.files(idx).default_roi_areas = handles.files(idx).roi_areas;
                handles.files(idx).default_roi_idcs = handles.files(idx).roi_idcs;
                handles.files(idx).loaded = true;
            catch
                warning('Import error - .mat file is not of the correct format, cannot load')
            end
        case {'.tif' '.tiff'}
            handles.files(idx).image = imread(handles.files(idx).filepath);
            handles.files(idx).dims = size(handles.files(idx).image);
            handles.files(idx).loaded = true;
            handles.files(idx).default_image = handles.files(idx).image;
    end
    handles.displayed_image = idx;
    initialize_displays(hObject,handles);
    update_displays(hObject,handles);
end

function handles = set_gui_defaults(hObject,handles)
handles.selected_rois = [];
handles.overlapping_rois = [nan nan];
handles.p_overlaps = [nan];
guidata(hObject,handles);

function select_rois(hObject,eventdata)
handles = guidata(hObject);
if eventdata.Button == 1 && ~isempty(handles.files(1).rois_ibar) && ~isempty(handles.files(2).rois_ibar)
    file_idcs = 1:numel(handles.files);
    idx = sub2ind(handles.files(1).dims,round(eventdata.IntersectionPoint(2)),round(eventdata.IntersectionPoint(1)));
    handles.selected_rois = [];
    for i = 1:numel(handles.files)
        if ~isempty(handles.files(i).rois_ibar(handles.files(i).rois_ibar(:,1) == idx,end))
            handles.selected_rois = [handles.selected_rois ; ...
                                     (i + (0*handles.files(i).rois_ibar(handles.files(i).rois_ibar(:,1) == idx,end))) ...
                                     handles.files(i).rois_ibar(handles.files(i).rois_ibar(:,1) == idx,end)];
        end
    end
    if ~isempty(handles.selected_rois)
        num_selected = size(handles.selected_rois,1);
        centroids_yx = zeros(num_selected,2);
        for r = 1:num_selected
            %centroids_yx(r,:) = handles.files(handles.selected_rois(r,1)).centroids(handles.selected_rois(r,2),:);
            % Added 20170511 HD to fix roi indexing
            centroids_yx(r,:) = handles.files(handles.selected_rois(r,1)).centroids(handles.files(handles.selected_rois(r,1)).roi_idcs == handles.selected_rois(r,2),:);
        end
        [~,min_dist_idx] = eudist(centroids_yx,eventdata.IntersectionPoint([2 1]));
        handles.selected_rois = handles.selected_rois(min_dist_idx,:);
        this_file = file_idcs(file_idcs==handles.selected_rois(1));
        other_file = file_idcs(file_idcs~=handles.selected_rois(1));
        [~,overlap_candidate] = eudist(handles.files(other_file).centroids,centroids_yx(min_dist_idx,:));
        
        overlap_candidate = handles.files(other_file).roi_idcs(overlap_candidate); % Added 20170511 HD to correct roi indexing
        
        these_idcs = handles.files(this_file).rois_ibar(:,4)==handles.selected_rois(2);
        other_idcs = handles.files(other_file).rois_ibar(:,4)==overlap_candidate;
        
        overlap_p = calculate_overlap(handles.files(this_file).rois_ibar(these_idcs,1),handles.files(other_file).rois_ibar(other_idcs,1));
        if overlap_p > 0  
            if (handles.overlapping_rois(:,other_file) == overlap_candidate) == (handles.overlapping_rois(:,this_file) == handles.selected_rois(2))
                if max(handles.overlapping_rois(:,other_file) == overlap_candidate) == 1
                    toremove = handles.overlapping_rois(:,this_file) == handles.selected_rois(2);
                    handles.overlapping_rois(toremove,:) = [];
                    handles.p_overlaps(toremove) = [];
                else
                    handles.overlapping_rois = [handles.overlapping_rois ; 0 0];
                    handles.overlapping_rois(end,[this_file other_file]) = [handles.selected_rois(2) overlap_candidate];
                    handles.p_overlaps = [handles.p_overlaps ; overlap_p];
                end
            else
                prev_selected = [find(handles.overlapping_rois(:,this_file) == handles.selected_rois(2)) ; ...
                                 find(handles.overlapping_rois(:,other_file) == overlap_candidate)];
                for i = 1:size(prev_selected,1)
                    for j = 1:size(handles.overlapping_rois,2)
                        prev_roi = handles.overlapping_rois(prev_selected(i),j);
                        prev_idcs = handles.files(j).rois_ibar(:,4)==prev_roi;
                        handles.files(j).roi_status(prev_idcs,:) = [1+(0*handles.files(j).roi_status(prev_idcs,1)) (0*handles.files(j).roi_status(prev_idcs,2))];
                    end
                end
                handles.overlapping_rois(prev_selected,:) = [];
                handles.overlapping_rois = [handles.overlapping_rois ; 0 0];
                handles.overlapping_rois(end,[this_file other_file]) = [handles.selected_rois(2) overlap_candidate];
                handles.p_overlaps(prev_selected) = [];
                handles.p_overlaps = [handles.p_overlaps ; overlap_p];
            end
            handles.files(this_file).roi_status(these_idcs,:) = ~handles.files(this_file).roi_status(these_idcs,:);
            handles.files(other_file).roi_status(other_idcs,:) = ~handles.files(other_file).roi_status(other_idcs,:);
            update_displays(hObject,handles);
            if handles.targets_overlapped
                handles = detect_target_overlap(handles);
            end
        end
    end
end
guidata(hObject,handles)

function handles = reload_roi_defaults(handles)
for i = 1:numel(handles.files)
    handles.files(i).image = handles.files(i).default_image;
    handles.files(i).dims = handles.files(i).dims;
    handles.files(i).rois_ibar = handles.files(i).default_rois_ibar;
    handles.files(i).centroids = handles.files(i).default_centroids;
    handles.files(i).roi_idcs = handles.files(i).default_roi_idcs;
    handles.files(i).roi_areas = handles.files(i).default_roi_areas;
    handles.files(i).num_rois = numel(handles.files(i).roi_idcs);
    handles.files(i).roi_status = [ones(size(handles.files(i).rois_ibar,1),1) zeros(size(handles.files(i).rois_ibar,1),1)];
end
handles.overlapping_rois = [nan nan];
handles.p_overlaps = [nan];

function handles = lasso_targets(hObject,~,~)
handles = guidata(hObject);
if size(handles.tformed_targets_yx,1) > 1
    if max(handles.selected_targets) == 0
        try
            handles.plotted_targets.CData(:,2) = 0;
            lasso_selected = imfreehand(gca);
            lasso_area = wait(lasso_selected);
            handles.selected_targets = 0 * handles.selected_targets;
            handles.selected_targets = inpolygon(handles.tformed_targets_yx(:,2),handles.tformed_targets_yx(:,1),lasso_area(:,1),lasso_area(:,2));
            delete(lasso_selected);
            handles.plotted_targets.CData(:,2) = handles.selected_targets;
        end
    else
        p = get(gca,'CurrentPoint');
        s = handles.selected_targets;
        yx = [handles.plotted_targets.YData' handles.plotted_targets.XData'];
        set(gcf,'WindowButtonMotionFcn',{@drag, handles, p, yx, s})
    end
    guidata(hObject,handles);
end

function drag(hObject,event,handles,p,yx,s)
np = get(gca,'CurrentPoint');
dyx = np(1,[2 1])- p(1,[2 1]);
nyx = bsxfun(@plus,yx,dyx);
% Added HD 20170510
nyx(nyx<=0) = 1;
nyx(nyx(:,1)>handles.files(handles.displayed_image).dims(1),1) = handles.files(handles.displayed_image).dims(1);
nyx(nyx(:,2)>handles.files(handles.displayed_image).dims(2),2) = handles.files(handles.displayed_image).dims(2);
handles.tformed_targets_yx(s,:) = nyx(s,:);
handles.plotted_targets.XData(s) = nyx(s,2);
handles.plotted_targets.YData(s) = nyx(s,1);
if handles.targets_overlapped
    handles = detect_target_overlap(handles);
end
guidata(hObject,handles);


function [handles] = plot_targets(hObject,handles)
for i = 1:numel(handles.dil_targets_handles)
    delete(handles.dil_targets_handles(i));
end
[handles.dil_targets_handles,handles.dil_targets_yx,~] = draw_rois_local(handles.tformed_targets_yx,handles.image_axis,handles.files(handles.displayed_image).dims,'Color',[0.6 0.6 0.6]);
handles.original_targets.XData = handles.tformed_targets_default(:,2);
handles.original_targets.YData = handles.tformed_targets_default(:,1);
handles.plotted_targets.XData = handles.tformed_targets_yx(:,2);
handles.plotted_targets.YData = handles.tformed_targets_yx(:,1);
handles.plotted_targets.SizeData = handles.target_sizes(handles.targets_s2p_rois(:,3)+1);
handles.plotted_targets.CData = repmat([1 0 0],size(handles.tformed_targets_yx,1),1);
for i = 1:numel(handles.targets)
    handles.targets(i).scatter.XData = handles.files(i).targets_overlap_yx(:,2);
    handles.targets(i).scatter.YData = handles.files(i).targets_overlap_yx(:,1);
end
guidata(hObject,handles);

function [handles] = reset_targets(hObject,handles);
for i = 1:numel(handles.dil_targets_handles)
    delete(handles.dil_targets_handles(i));
end
handles.dil_targets_handles = [];
handles.original_targets.XData = [];
handles.original_targets.YData = [];
handles.plotted_targets.XData = [];
handles.plotted_targets.YData = [];
handles.plotted_targets.SizeData = handles.target_sizes(1);
handles.plotted_targets.CData = [1 0 1];
for i = 1:numel(handles.targets)
    handles.targets(i).scatter.XData = [];
    handles.targets(i).scatter.YData = [];
end
guidata(hObject,handles);


function [handles] = transform_rois(handles,idx)
tformed_ibar = [];
handles.files(idx).centroids = [];
handles.files(idx).roi_areas = [];
for r = 1:handles.files(idx).num_rois
    blank_im = 0 * handles.files(idx).image{1}(:,:,1);
    r_idcs = handles.files(idx).rois_ibar(:,4) == handles.files(idx).roi_idcs(r) & handles.files(idx).rois_ibar(:,3) == 1; % "& handles.files(idx).rois_ibar(:,3) == 1" Added 20170517 HD to fix ROI growing post transform
    im_idcs = handles.files(idx).rois_ibar(r_idcs,1);
    blank_im(im_idcs) = 1;
    blank_im = imwarp(blank_im,handles.tform,'OutputView',imref2d(size(blank_im)));
    blank_im(blank_im>0) = 1;
    these_idcs = find(blank_im==1);
    % If ROI hasn't been translated out of the image
    if numel(these_idcs) > 0 % Added 20170511 HD to fix ROI indexing
        b = bwperim(blank_im);
        b = find(b);
        border = 0 * these_idcs;
        border(ismember(these_idcs,b)) = 1;
        area = ~border;      
        tformed_ibar = [tformed_ibar ; these_idcs border area handles.files(idx).roi_idcs(r) * ones(numel(these_idcs),1)];
        [y,x] = ind2sub(size(blank_im),tformed_ibar(end-numel(these_idcs)+1:end,1));
        handles.files(idx).centroids = [handles.files(idx).centroids ; mean([y x],1)];
        handles.files(idx).roi_areas = [handles.files(idx).roi_areas ; numel(these_idcs)];
    end
end
handles.files(idx).rois_ibar = tformed_ibar;
handles.files(idx).roi_status = [ones(size(handles.files(idx).rois_ibar(:,1),1),1) zeros(size(handles.files(idx).rois_ibar(:,1),1),1)];
handles.files(idx).roi_idcs = unique(handles.files(idx).rois_ibar(:,4),'stable');
handles.files(idx).num_rois = numel(handles.files(idx).roi_idcs);


function [rois_ibar,centroids,areas,s2p_iscell] = load_S2P_rois(dat,dims)
s2p_iscell = find([dat.stat(:).iscell]);
rois_ibar = zeros(sum(cellfun(@(x) numel(x),{dat.stat(s2p_iscell).ypix}')),4);
prev_id = 1;
centroids = zeros(numel(s2p_iscell),2);
areas = centroids(:,1);
%min_x = dat.ops.xrange(1);
%min_y = dat.ops.yrange(1);
for r = 1:numel(s2p_iscell)
    %centroids(r,:) = ROI_text1([dat.stat(s2p_iscell(r)).ypix+min_y dat.stat(s2p_iscell(r)).xpix+min_x],1);
    %these_idcs = sub2ind(dims,dat.stat(s2p_iscell(r)).ypix+min_y,dat.stat(s2p_iscell(r)).xpix+min_x);
    centroids(r,:) = mean([dat.ops.yrange(dat.stat(s2p_iscell(r)).ypix)' dat.ops.xrange(dat.stat(s2p_iscell(r)).xpix)'],1); % Added 20170511 HD to fix ROI indexing
    these_idcs = sub2ind(dims,dat.ops.yrange(dat.stat(s2p_iscell(r)).ypix),dat.ops.xrange(dat.stat(s2p_iscell(r)).xpix))'; % Added 20170511 HD to fix ROI indexing
    areas(r) = numel(these_idcs);
    blank_im = zeros(dims);
    blank_im(these_idcs) = 1;
    b = bwperim(blank_im);
    b = find(b);
    border = 0 * these_idcs;
    border(ismember(these_idcs,b)) = 1;
    area = ~border; %1 + (0 * these_idcs);
    %rois_ibar(prev_id:prev_id+numel(these_idcs)-1,:) = [these_idcs border area  (r + (0 * these_idcs))];
    rois_ibar(prev_id:prev_id+numel(these_idcs)-1,:) = [these_idcs border area  (s2p_iscell(r) + (0 * these_idcs))]; % Added 20170511 HD to fix ROI indexing
    prev_id = prev_id + numel(these_idcs);
end

function initialize_displays(hObject,handles)
handles.image_display.CData = 0 * handles.files(handles.displayed_image).image{1}(:,:,1);
handles.roi_display.CData = 0 * handles.files(handles.displayed_image).image{1}(:,:,1);
guidata(hObject,handles)

function update_displays(hObject,handles)
non_displayed = find([1:numel(handles.files)] ~= handles.displayed_image);
handles.roi_display.CData = 0 * handles.roi_display.CData;
for i = non_displayed
    if ~isempty(handles.files(i).rois_ibar)
        handles.roi_display.CData(handles.files(i).rois_ibar(handles.files(i).roi_status(:,1)==1,1)) = handles.files(i).rois_ibar(handles.files(i).roi_status(:,1)==1,2) * 20;
        handles.roi_display.CData(handles.files(i).rois_ibar(handles.files(i).roi_status(:,2)==1,1)) = handles.files(i).rois_ibar(handles.files(i).roi_status(:,2)==1,3) * 20;
    end
end
i = handles.displayed_image;
if ~isempty(handles.files(i).rois_ibar)
    handles.roi_display.CData(handles.files(i).rois_ibar(handles.files(i).roi_status(:,1)==1,1)) = ...
        handles.roi_display.CData(handles.files(i).rois_ibar(handles.files(i).roi_status(:,1)==1,1)) + handles.files(i).rois_ibar(handles.files(i).roi_status(:,1)==1,2) * 50;
    
    handles.roi_display.CData(handles.files(i).rois_ibar(handles.files(i).roi_status(:,1)==1,1)) = ...
        handles.roi_display.CData(handles.files(i).rois_ibar(handles.files(i).roi_status(:,1)==1,1)) + handles.files(i).rois_ibar(handles.files(i).roi_status(:,1)==1,2) * 50;
    
    handles.roi_display.CData(handles.files(i).rois_ibar(handles.files(i).roi_status(:,2)==1,1)) = ...
        handles.roi_display.CData(handles.files(i).rois_ibar(handles.files(i).roi_status(:,2)==1,1)) + handles.files(i).rois_ibar(handles.files(i).roi_status(:,2)==1,3) * 50;
end
if ~isempty(handles.files(handles.displayed_image).image)
    if handles.im_flag
        handles.image_display.CData = handles.files(handles.displayed_image).image{handles.proc_flag}(:,:,handles.displayed_mask);
        lims = [min(handles.image_display.CData(:)) max(handles.image_display.CData(:))];
        if lims(1) ~= lims(2)
            handles.image_display.Parent.CLim = [min(handles.image_display.CData(:)) max(handles.image_display.CData(:))];
        else
            handles.image_display.Parent.CLim = [lims(1)-1 lims(1)+1];
        end
    else
        if ~isempty(handles.files(i).rois_ibar)
            handles.image_display.Parent.CLim = [0 10];
            handles.image_display.CData = 0 * handles.image_display.CData;
            handles.image_display.CData(handles.files(i).rois_ibar(handles.files(i).roi_status(:,1)==1,1)) = ...
                handles.image_display.CData(handles.files(i).rois_ibar(handles.files(i).roi_status(:,1)==1,1)) + handles.files(i).rois_ibar(handles.files(i).roi_status(:,1)==1,2) * 8;
            
            handles.image_display.CData(handles.files(i).rois_ibar(handles.files(i).roi_status(:,2)==1,1)) = ...
                handles.image_display.CData(handles.files(i).rois_ibar(handles.files(i).roi_status(:,2)==1,1)) + handles.files(i).rois_ibar(handles.files(i).roi_status(:,2)==1,3) * 8;
            handles.image_display.CData(handles.image_display.CData>10) = 10;
        end
    end
end
guidata(hObject,handles);

function [distances,min_dist_idx] = eudist(points2search,point2match)
distances = sqrt(sum((points2search - repmat(point2match,size(points2search,1),1)) .^2,2));
[~,order] = sort(distances,'Ascend');
min_dist_idx = order(1);

function p_overlap = calculate_overlap(roi_idcs1,roi_idcs2)
%p_overlap = numel(intersect(roi_idcs1,roi_idcs2)) / numel(union(roi_idcs1,roi_idcs2));
p_overlap = numel(intersect(roi_idcs1,roi_idcs2)) / numel(roi_idcs1);

function [overlap_idcs,overlaps] = overlap_detector(target_centroids,roi_centroids,overlap_threshold)
overlap_idcs = nan(size(target_centroids,1),1);
overlaps = nan(size(target_centroids,1),1);
for c = 1:size(target_centroids,1)
    [d,mi] = eudist(roi_centroids,target_centroids(c,:));
    overlaps(c) = d(mi);
    if overlaps(c) <= overlap_threshold
        present_flag = overlap_idcs == mi;
        if max(present_flag) == 0
            overlap_idcs(c) = mi;
        else
            if overlaps(present_flag) > overlaps(c)
                overlap_idcs(present_flag) = nan;
                overlap_idcs(c) = mi;
            end
        end
    end
end

function handles = detect_target_overlap(handles)
if handles.targets_loaded && (handles.files(1).loaded || handles.files(2).loaded)
    handles.targets_s2p_rois = nan * handles.targets_s2p_rois;
    [idcs,~] = overlap_detector(handles.tformed_targets_yx,handles.files(1).centroids,handles.targets_overlap_thold);
    handles.targets_s2p_rois(~isnan(idcs),1) = handles.files(1).roi_idcs(idcs(~isnan(idcs)));
    for i = 1:numel(handles.files)
        handles.targets(i).scatter.XData = [];
        handles.targets(i).scatter.YData = [];
    end
    for i = 1:size(handles.targets_s2p_rois,1)
        if ~isnan(handles.targets_s2p_rois(i,1))
            id = handles.files(1).rois_ibar(:,end)==handles.targets_s2p_rois(i,1) & handles.files(1).rois_ibar(:,2)==1;
            [y,x] = ind2sub(handles.files(1).dims,handles.files(1).rois_ibar(id,1));
            handles.targets(1).scatter.XData = [handles.targets(1).scatter.XData x'];
            handles.targets(1).scatter.YData = [handles.targets(1).scatter.YData y'];
            overlapping_roi = handles.overlapping_rois(handles.overlapping_rois(:,1)==handles.targets_s2p_rois(i,1),2);
            if ~isempty(overlapping_roi)
                handles.targets_s2p_rois(i,2) = overlapping_roi;
                id = handles.files(2).rois_ibar(:,end)==overlapping_roi & handles.files(2).rois_ibar(:,2)==1;
                [y,x] = ind2sub(handles.files(2).dims,handles.files(2).rois_ibar(id,1));
                handles.targets(2).scatter.XData = [handles.targets(2).scatter.XData x'];
                handles.targets(2).scatter.YData = [handles.targets(2).scatter.YData y'];
            end
        end
    end
    for i = 1:numel(handles.files)
        handles.files(i).targets_overlap_yx = [handles.targets(i).scatter.YData' handles.targets(i).scatter.XData'];
    end
    handles.targets_s2p_rois(:,3) = ~isnan(sum(handles.targets_s2p_rois(:,1:2),2));
    handles.plotted_targets.SizeData = handles.target_sizes(handles.targets_s2p_rois(:,3)+1);
    handles.targets_overlapped = true; 
end


function handles = reset_targets_overlap(handles)
handles.targets_s2p_rois(:,1:2) = nan * handles.targets_s2p_rois(:,1:2);
handles.targets_s2p_rois(:,3) = 0 * handles.targets_s2p_rois(:,3);
handles.plotted_targets.SizeData = handles.target_sizes(1);
for i = 1:numel(handles.files)
    handles.files(i).targets_overlap_yx = [nan nan];
    handles.targets(i).scatter.XData = [];
    handles.targets(i).scatter.YData = [];
end
handles.targets_overlapped = false;

function handles = reset_rois_overlap(handles)
for idx = 1:numel(handles.files)
    handles.files(idx).roi_status = [ones(size(handles.files(idx).rois_ibar(:,1),1),1) zeros(size(handles.files(idx).rois_ibar(:,1),1),1)];
end

function [handles,outlines,all_mask] = draw_rois_local(rois,axis_handle,dims,varargin);
r = 10;
col = [0.7 0.7 0.7];
wid = 1;
plot_flag = 1;
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'radius')
       r = varargin{v+1};
    elseif strcmpi(varargin{v},'Color')
       col = varargin{v+1};
    elseif strcmpi(varargin{v},'LineWidth')
        wid = varargin{v+1};
    elseif strcmpi(varargin{v},'Plot')
        plot_flag = varargin{v+1}; 
    end
end
num_conts = size(rois,1);
if size(rois,2) == 1
   temp = rois;
   rois = [];
   [rois(:,1),rois(:,2)] = ind2sub([dims(1) dims(2)],temp);
end
x_mask = zeros(dims(1),dims(2),num_conts);
blank_im = zeros(dims(1),dims(2));
th = 0:pi/50:2*pi;
for j = 1:num_conts
    yxunit(:,2) = round(r * cos(th) + rois(j,2));
    yxunit(:,1) = round(r * sin(th) + rois(j,1));
    yxunit(yxunit <= 0) = 1;
    yxunit(yxunit(:,1)>dims(1),1) = dims(1);
    yxunit(yxunit(:,2)>dims(2),2) = dims(2);
    z = blank_im; z(sub2ind(size(z),yxunit(:,1),yxunit(:,2))) = 1;
    x_mask(:,:,j) = imfill(z);
end
all_mask = sum(x_mask,3); 
all_mask(all_mask>0) = 1;
[B,~] = bwboundaries(all_mask);
num_outlines = numel(B);
outlines = cell(num_outlines,1);
handles = zeros(num_outlines,1);
for i = 1:num_outlines
    x = downsample(B{i}(:,2),3);
    x = [x(:) ; x(1)];
    y = downsample(B{i}(:,1),3);
    y = [y(:) ; y(1)];
    if plot_flag
        hold on
        handles(i) = plot(axis_handle,x,y,'Color',col,'LineWidth',wid);
    end
    outlines{i} = [x' y'];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ButtonDownFcns %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on key press with focus on load_fixed_button and none of its controls.
function load_fixed_button_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to load_fixed_button (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

figure1_KeyPressFcn(hObject, eventdata, handles)


% --- Executes on key press with focus on load_transform_button and none of its controls.
function load_transform_button_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to load_transform_button (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
figure1_KeyPressFcn(hObject, eventdata, handles)


% --- Executes on key press with focus on load_moving_button and none of its controls.
function load_moving_button_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to load_moving_button (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
figure1_KeyPressFcn(hObject, eventdata, handles)


% --- Executes on key press with focus on load_targets_button and none of its controls.
function load_targets_button_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to load_targets_button (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
figure1_KeyPressFcn(hObject, eventdata, handles)


% --- Executes on key press with focus on save_tform_button and none of its controls.
function save_tform_button_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to save_tform_button (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
figure1_KeyPressFcn(hObject, eventdata, handles)


% --- Executes on key press with focus on register_button and none of its controls.
function register_button_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to register_button (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
figure1_KeyPressFcn(hObject, eventdata, handles)


% --- Executes on key press with focus on reset_button and none of its controls.
function reset_button_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to reset_button (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
figure1_KeyPressFcn(hObject, eventdata, handles)

% --- Executes on key press with focus on save_button and none of its controls.
function save_button_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
figure1_KeyPressFcn(hObject, eventdata, handles)


% --- Executes on button press in detect_overlap_button.
function detect_overlap_button_Callback(hObject, eventdata, handles)
% hObject    handle to detect_overlap_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.files(1).rois_ibar) && ~isempty(handles.files(2).rois_ibar)
    for i = 1:numel(handles.files)
        handles.files(i).roi_status = [1+(0*handles.files(i).roi_status(:,1)) (0*handles.files(i).roi_status(:,2))];
    end
    handles.overlapping_rois = [nan nan];
    handles.p_overlaps = nan;
    update_displays(hObject,handles)
    overlap_im = zeros(handles.files(1).dims);
    for i = 1:numel(handles.files)
        overlap_im(unique(handles.files(i).rois_ibar(:,1))) = overlap_im(unique(handles.files(i).rois_ibar(:,1))) + 1;
    end
    overlap_idcs = find(overlap_im == 2);
    overlapping_rois1 = unique(handles.files(1).rois_ibar(ismember(handles.files(1).rois_ibar(:,1),overlap_idcs),4));
    for j = 1:numel(overlapping_rois1)
        roi_idcs1 = handles.files(1).rois_ibar(handles.files(1).rois_ibar(:,4)==overlapping_rois1(j),1);
        overlapping_rois2 = unique(handles.files(2).rois_ibar(ismember(handles.files(2).rois_ibar(:,1),roi_idcs1),4));
        roi_idcs2 = cell(numel(overlapping_rois2),1);
        p_overlap = zeros(numel(overlapping_rois2),1);
        for i = 1:numel(overlapping_rois2)
            roi_idcs2{i} = handles.files(2).rois_ibar(handles.files(2).rois_ibar(:,4)==overlapping_rois2(i),1);
            p_overlap(i) = calculate_overlap(roi_idcs1,roi_idcs2{i});
        end
        [~,max_p] = max(p_overlap);
        p_overlap = p_overlap(max_p);
        overlapping_rois2 = overlapping_rois2(max_p);
        if p_overlap > handles.overlap_thold
            if max(handles.overlapping_rois(:,2) == overlapping_rois2) ~= 1
                handles.overlapping_rois = [handles.overlapping_rois ; overlapping_rois1(j) overlapping_rois2];
                handles.p_overlaps = [handles.p_overlaps ; p_overlap];
                for i = 1:numel(handles.files)
                    handles.files(i).roi_status(handles.files(i).rois_ibar(:,4)==handles.overlapping_rois(end,i),:) = ~handles.files(i).roi_status(handles.files(i).rois_ibar(:,4)==handles.overlapping_rois(end,i),:);
                end
            else
                if p_overlap > handles.p_overlaps(handles.overlapping_rois(:,2) == overlapping_rois2)
                    overlap2remove = handles.overlapping_rois(handles.overlapping_rois(:,2) == overlapping_rois2,:);
                    for i = 1:numel(handles.files)
                        handles.files(i).roi_status(handles.files(i).rois_ibar(:,4)==overlap2remove(i),:) = ~handles.files(i).roi_status(handles.files(i).rois_ibar(:,4)==overlap2remove(i),:);
                    end
                    remove_flag = handles.overlapping_rois(:,2) == overlapping_rois2;
                    handles.overlapping_rois(remove_flag,:) = []; % Added HD 20170502
                    handles.p_overlaps(remove_flag) = [];
                    handles.overlapping_rois = [handles.overlapping_rois ; overlapping_rois1(j) overlapping_rois2];
                    handles.p_overlaps = [handles.p_overlaps ; p_overlap];
                    for i = 1:numel(handles.files)
                        handles.files(i).roi_status(handles.files(i).rois_ibar(:,4)==handles.overlapping_rois(end,i),:) = ~handles.files(i).roi_status(handles.files(i).rois_ibar(:,4)==handles.overlapping_rois(end,i),:);
                    end
                end
            end
        end
    end
    if handles.targets_overlapped
        handles = detect_target_overlap(handles);
    end
    update_displays(hObject,handles);
end
guidata(hObject,handles);

% --- Executes on button press in detect_target_overlap_button.
function detect_target_overlap_button_Callback(hObject, eventdata, handles)
% hObject    handle to detect_target_overlap_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = detect_target_overlap(handles);
guidata(hObject,handles);



% --- Executes on button press in reset_targets_overlap_button.
function reset_targets_overlap_button_Callback(hObject, eventdata, handles)
% hObject    handle to reset_targets_overlap_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles.targets_s2p_rois = nan * handles.targets_s2p_rois;
% handles.plotted_targets.SizeData = handles.target_sizes(1);
% for i = 1:numel(handles.files)
%     handles.targets(i).scatter.XData = [];
%     handles.targets(i).scatter.YData = [];
% end
% handles.targets_overlapped = false;
handles = reset_targets_overlap(handles);
guidata(hObject,handles);

function overlap_threshold_input_Callback(hObject, eventdata, handles)
% hObject    handle to overlap_threshold_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of overlap_threshold_input as text
%        str2double(get(hObject,'String')) returns contents of overlap_threshold_input as a double
thold = str2double(get(hObject,'String'));
if thold > 1
    thold = 1;
    set(hObject,'String',num2str(thold));
elseif thold <= 0
    set(hObject,'String',handles.overlap_thold)
elseif isnan(thold)
    set(hObject,'String',handles.overlap_thold)
end
handles.overlap_thold = str2double(get(hObject,'String'));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function overlap_threshold_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to overlap_threshold_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in reset_overlap_button.
function reset_overlap_button_Callback(hObject, eventdata, handles)
% hObject    handle to reset_overlap_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.files(1).rois_ibar) && ~isempty(handles.files(2).rois_ibar)
    for i = 1:numel(handles.files)
        handles.files(i).roi_status = [1+(0*handles.files(i).roi_status(:,1)) (0*handles.files(i).roi_status(:,2))];
    end
    handles.overlapping_rois = [nan nan];
    handles.p_overlaps = nan;
    update_displays(hObject,handles)
    if handles.targets_overlapped
        handles = detect_target_overlap(handles);
    end
end
guidata(hObject,handles)


% --- Executes on key press with focus on detect_overlap_button and none of its controls.
function detect_overlap_button_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to detect_overlap_button (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
figure1_KeyPressFcn(hObject, eventdata, handles)


% --- Executes on key press with focus on reset_overlap_button and none of its controls.
function reset_overlap_button_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to reset_overlap_button (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
figure1_KeyPressFcn(hObject, eventdata, handles)


% --- Executes on key press with focus on overlap_threshold_input and none of its controls.
function overlap_threshold_input_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to overlap_threshold_input (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
%figure1_KeyPressFcn(hObject, eventdata, handles)


% --- Executes on button press in save_session_button.
function save_session_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_session_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.files(1).filepath) || ~isempty(handles.files(2).filepath)
    name_idcs = find(cellfun(@(x) ~isempty(x),{handles.files(:).filepath}));
    name_chunks = cell(1,numel(name_idcs));
    ff = name_chunks;
    for i = 1:numel(name_idcs)
        [dd,ff{i}] = fileparts(handles.files(name_idcs(i)).filepath);
        name_chunks{i} = ff{i}(3:end);
        name_chunks{i} = strrep(name_chunks{i},'proc','');
    end
    default_name = [name_chunks{:}];
    default_name(end) = [];
    default_name = [default_name '_session'];
    [f,d] = uiputfile([dd filesep '*.mat'],'Save...',[dd filesep default_name '.mat']);
    if d ~= 0
        fnames = fieldnames(handles);
        s_idx = find(cellfun(@(x) strcmp(x,'files'),fnames));
        out_handles = [];
        for fn = s_idx:numel(fnames)
            out_handles.(fnames{fn}) = handles.(fnames{fn});
        end
        save([d filesep f],'out_handles');
    end
else
    uiwait(msgbox('No data loaded - please load at least 1 session before trying to save', 'Error','error'));
end


% --- Executes on button press in load_session_button.
function load_session_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_session_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[f,d] = uigetfile('*.mat');
if d ~= 0
    load([d filesep f])
    fnames = fieldnames(handles);
    in_fnames = fieldnames(out_handles);
    if min(isfield(handles,in_fnames)) == 1
        s_idx = find(cellfun(@(x) strcmp(x,'files'),fnames));
        for fn = s_idx:numel(fnames)
            handles.(fnames{fn}) = out_handles.(fnames{fn});
        end
        if ~isempty(handles.files(1).filepath)
            handles.fixed_textbox.String = handles.files(1).filepath;
        else
            handles.fixed_textbox.String = '...';
        end
        if ~isempty(handles.files(2).filepath)
            handles.moving_textbox.String = handles.files(2).filepath;
        else
            handles.moving_textbox.String = '...';
        end
        if ~isempty(handles.targets_filepath)
            handles.targets_textbox.String = handles.targets_filepath;
        else
            handles.targets_textbox.String = '...';
        end
        if ~isempty(handles.transform_filepath)
            handles.transform_textbox.String = handles.transform_filepath;
        else
            handles.transform_textbox.String = '...';
        end
        handles.overlap_threshold_input.String = handles.overlap_thold;
    else
        sprintf(['Error loading session'])
    end
    [handles] = reset_targets(hObject,handles);
    initialize_displays(hObject,handles)
    update_displays(hObject,handles);
    if handles.targets_loaded
        handles = plot_targets(hObject,handles);
    end
    guidata(hObject,handles);
end


function target_overlap_threshold_input_Callback(hObject, eventdata, handles)
% hObject    handle to target_overlap_threshold_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of target_overlap_threshold_input as text
%        str2double(get(hObject,'String')) returns contents of target_overlap_threshold_input as a double
thold = str2double(get(hObject,'String'));
if isnan(thold)
    set(hObject,'String',handles.targets_overlap_thold)
end
handles.targets_overlap_thold = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function target_overlap_threshold_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to target_overlap_threshold_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on key press with focus on save_session_button and none of its controls.
function save_session_button_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to save_session_button (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
figure1_KeyPressFcn(hObject, eventdata, handles)


% --- Executes on key press with focus on load_session_button and none of its controls.
function load_session_button_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to load_session_button (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
figure1_KeyPressFcn(hObject, eventdata, handles)


% --- Executes on key press with focus on detect_target_overlap_button and none of its controls.
function detect_target_overlap_button_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to detect_target_overlap_button (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
figure1_KeyPressFcn(hObject, eventdata, handles)


% --- Executes on key press with focus on reset_targets_overlap_button and none of its controls.
function reset_targets_overlap_button_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to reset_targets_overlap_button (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
figure1_KeyPressFcn(hObject, eventdata, handles)


% --- Executes on key press with focus on target_overlap_threshold_input and none of its controls.
function target_overlap_threshold_input_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to target_overlap_threshold_input (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
%figure1_KeyPressFcn(hObject, eventdata, handles)


% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
if handles.overlay
    order = [1 0];
    C = get (handles.image_axis, 'CurrentPoint');
    C = C(1,1:2);
    if (C(1) > handles.files(1).dims(1) | C(2) > handles.files(1).dims(2)) || min(C)<0
        handles.cursors(2).scatter.XData = nan;
        handles.cursors(2).scatter.YData = nan;
        C = get (handles.roi_axis, 'CurrentPoint');
        C = C(1,1:2);
        idcs = (~order)+1;
    else
        idcs = order+1;
    end
    if (C(1) <= handles.files(idcs(1)).dims(1) | C(2) <= handles.files(idcs(1)).dims(2)) || min(C)>0
        handles.cursors(idcs(2)).scatter.XData = nan;
        handles.cursors(idcs(2)).scatter.YData = nan;
        handles.cursors(idcs(1)).scatter.XData = C(1,1);
        handles.cursors(idcs(1)).scatter.YData = C(1,2);
    end
end


% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if max(handles.selected_targets) == 1
    set(gcf,'WindowButtonMotionFcn',@figure1_WindowButtonMotionFcn)
end

% --- Executes on mouse press over figure background.
function figure1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text1.
function text1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to text1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over target_overlap_threshold_input.
function target_overlap_threshold_input_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to target_overlap_threshold_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
