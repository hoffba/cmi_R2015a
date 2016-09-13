function varargout = setRegions(varargin)
% SETREGIONS MATLAB code for setRegions.fig
%      SETREGIONS, by itself, creates a new SETREGIONS or raises the existing
%      singleton*.
%
%      H = SETREGIONS returns the handle to a new SETREGIONS or the handle to
%      the existing singleton*.
%
%      SETREGIONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SETREGIONS.M with the given input arguments.
%
%      SETREGIONS('Property','Value',...) creates a new SETREGIONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before setRegions_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to setRegions_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help setRegions

% Last Modified by GUIDE v2.5 01-Aug-2016 10:04:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @setRegions_OpeningFcn, ...
                   'gui_OutputFcn',  @setRegions_OutputFcn, ...
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
% --- Executes just before setRegions is made visible.
function setRegions_OpeningFcn(hObject, ~, handles, varargin)
handles.output = hObject;
handles.i = 1;
handles.thi = [];
handles.cuti = [];
handles.normchk = true;
handles.Ci = [];
handles.maxnp = [];
handles.filtchk = true;
handles.filttype = 'Wiener';
handles.filtstr = [3,3];
handles.SPopts = struct('Xvec',{nan},'Yvec',{nan},'Xmin',{nan},'Ymin',{nan},...
                        'Xmax',{nan},'Ymax',{nan});

% Find saved default PRM settings
fpath = fullfile(fileparts(which('cmi')),'PRMdefs');
fnames = dir(fullfile(fpath,'PRMdef*.mat'));
handles.defs = struct('name','','thresh',[],'cutoff',[],'cmap',[],...
                      'prmmap',[],'SPopts',[]);
chk = false(1,length(fnames));
for i = 1:length(fnames)
    tstr = load(fullfile(fpath,fnames(i).name));
    if all(isfield(tstr,{'name','thresh','cutoff','cmap','prmmap'}))
        if ~isfield(tstr,'SPopts')
            tstr.SPopts = handles.SPopts;
        end
        handles.defs(i) = orderfields(tstr,handles.defs);
    else
        chk(i) = true;
    end
end
handles.defs(chk) = [];
set(handles.popup_defs,'String',{'Select Default Parameters:',handles.defs(:).name});

guidata(hObject, handles);
% --- Outputs from this function are returned to the command line.
function varargout = setRegions_OutputFcn(~, ~, handles) 
varargout{1} = handles.output;

% --- Executes on button press in button_color.
function button_color_Callback(~, ~, handles) %#ok<*DEFNU>
handles.cmap(handles.i,:) = uisetcolor(handles.cmap(handles.i,:));
set(handles.axes_color,'Color',handles.cmap(handles.i,:));
guidata(gcbo,handles);

% --- Executes on selection change in list_regions.
function list_regions_Callback(~, ~, handles)
i = get(handles.list_regions,'Value');
if isempty(i)
    set(handles.list_regions,'Value',handles.i);
else
    i = i(1);
    handles.i = i;
    guidata(gcbo,handles);
    set(handles.edit_name,'String',handles.prmmap{i,2});
    set(handles.table_C,'Data',num2cell(handles.prmmap{handles.i,1})');
    set(handles.axes_color,'Color',handles.cmap(i,:));
end

% --- Changes current region's name
function edit_name_Callback(~,~,handles)
str = get(handles.edit_name,'String');
if ~isempty(str)
    handles.prmmap{handles.i,2} = str;
    set(handles.list_regions,'String',handles.prmmap(:,2));
    guidata(gcbo,handles);
else
    set(handles.edit_name,'String',handles.prmmap{handles.i,2});
end

% --- Executes when entered data in editable cell(s) in table_C.
function table_C_CellEditCallback(~, ~, handles)
handles.prmmap{handles.i,1} = cell2mat(get(handles.table_C,'Data'))';
guidata(gcbo,handles);

% --- Combines selected regions.
function button_combine_Callback(~, ~, handles)
i = get(handles.list_regions,'Value');
if length(i)>1
    handles.prmmap{i(1),1} = cell2mat(handles.prmmap(i,1));
    handles.prmmap(i(2:end),:) = [];
    handles.cmap(i(2:end),:) = [];
    set(handles.list_regions,'Value',handles.i,'String',handles.prmmap(:,2))
    set(handles.table_C,'Data',num2cell(handles.prmmap{i(1),1})');
    guidata(gcbo,handles);
end

% --- Executes on button press in button_new.
function button_new_Callback(~, ~, handles)
nr = size(handles.prmmap,1) + 1;
handles.prmmap(nr,:) = {false(1,size(handles.thresh,1)),'NewRegion'};
handles.i = nr;
handles.cmap(nr,:) = [0,0,0];
set(handles.list_regions,'String',handles.prmmap(:,2),'Value',nr);
set(handles.axes_color,'Color',handles.cmap(nr,:));
set(handles.table_C,'Data',num2cell(handles.prmmap{nr,1})');
set(handles.edit_name,'String',handles.prmmap{nr,2});
guidata(gcbo,handles);

% --- Executes on button press in button_remove.
function button_remove_Callback(~, ~, handles)
i = get(handles.list_regions,'Value');
handles.cmap(i,:) = [];
handles.prmmap(i,:) = [];
handles.i = 1;
set(handles.list_regions,'Value',handles.i,'String',handles.prmmap(:,2));
set(handles.edit_name,'String',handles.prmmap(handles.i,2));
set(handles.table_C,'Data',num2cell(handles.prmmap{handles.i,1})');
set(handles.axes_color,'Color',handles.cmap(handles.i,:));
guidata(gcbo,handles);

% --- Saves region options.
function button_ok_Callback(~, ~, ~)
uiresume(gcbf);

% --- Executes on button press in button_cancel.
function button_cancel_Callback(~, ~, handles)
close(handles.figure1);

% --- Select default PRM options from popup menu.
function popup_defs_Callback(~, ~, handles)
n = get(handles.popup_defs,'Value');
if n>1
    
    handles.thresh = handles.defs(n-1).thresh;
    handles.cutoff = handles.defs(n-1).cutoff;
    handles.cmap = handles.defs(n-1).cmap;
    handles.prmmap = handles.defs(n-1).prmmap;
    handles.SPopts = handles.defs(n-1).SPopts;
    handles.i = 1;
    
    set(handles.table_thresh,'Data',num2cell(handles.thresh));
    set(handles.table_crop,'Data',num2cell(handles.cutoff));
    set(handles.table_C,'Data',num2cell(handles.prmmap{handles.i,1})');
    set(handles.list_regions,'Value',handles.i,'String',handles.prmmap(:,2));
    set(handles.edit_name,'String',handles.prmmap{handles.i,2});
    set(handles.axes_color,'Color',handles.cmap(handles.i,:));
    
    set(handles.edit_SPdimX,'String',handles.SPopts.Xvec);
    set(handles.edit_SPdimY,'String',handles.SPopts.Yvec);
    set(handles.edit_SPminX,'String',handles.SPopts.Xmin);
    set(handles.edit_SPmaxX,'String',handles.SPopts.Xmax);
    set(handles.edit_SPminY,'String',handles.SPopts.Ymin);
    set(handles.edit_SPmaxY,'String',handles.SPopts.Ymax);
    set(handles.edit_maxnp,'String',handles.SPopts.Nmax);
    set(handles.checkbox_showScatter,'Value',handles.SPopts.show);
    
    guidata(gcbo,handles);
end

% --- Executes when selected cell(s) is changed in table_crop.
function table_crop_CellSelectionCallback(~, edata, handles)
if ~isempty(edata.Indices)
    handles.cuti = edata.Indices(1);
else
    handles.cuti = [];
end
guidata(gcbo,handles);

% --- Executes when selected cell(s) is changed in table_thresh.
function table_thresh_CellSelectionCallback(~, edata, handles)
if ~isempty(edata.Indices)
    handles.thi = edata.Indices(1);
else
    handles.thi = [];
end
guidata(gcbo,handles);

% --- Executes when selected cell(s) is changed in table_C.
function table_C_CellSelectionCallback(~, edata, handles)
if ~isempty(edata.Indices)
    handles.Ci = edata.Indices(2); % select Col
else
    handles.Ci = [];
end
guidata(gcbo,handles);

% --- Executes when entered data in editable cell(s) in table_thresh.
function table_thresh_CellEditCallback(~, ~, handles)
handles.thresh = cell2mat(get(handles.table_thresh,'Data'));
tdims = handles.thresh(:,1:2);
mxdim = size(get(handles.list_imgs,'String'),1);
tdims(tdims<0) = 0;
tdims(tdims>mxdim) = mxdim;
handles.thresh(:,1:2) = tdims;
guidata(gcbo,handles);

% --- Executes when entered data in editable cell(s) in table_crop.
function table_crop_CellEditCallback(~, ~, handles)
handles.cutoff = cell2mat(get(handles.table_crop,'Data'));
tdims = handles.cutoff(:,1);
mxdim = size(get(handles.list_imgs,'String'),1);
tdims(tdims<0) = 0;
tdims(tdims>mxdim) = mxdim;
handles.cutoff(:,1) = tdims;
guidata(gcbo,handles);

% --- Executes on button press in button_th_add.
function button_th_add_Callback(~, ~, handles)
handles.thresh(end+1,:) = handles.thresh(end,:);
for i = 1:size(handles.prmmap,1)
    handles.prmmap{i,1}(:,end+1) = handles.prmmap{i,1}(:,end);
end
set(handles.table_C,'Data',num2cell(handles.prmmap{handles.i,1})');
set(handles.table_thresh,'Data',num2cell(handles.thresh));
guidata(gcbo,handles);

% --- Executes on button press in button_th_del.
function button_th_del_Callback(~, ~, handles)
if ~isempty(handles.thi);
    handles.thresh(handles.thi,:) = [];
    for i = 1:size(handles.prmmap,1)
        handles.prmmap{i,1}(:,end) = [];
    end
    set(handles.table_C,'Data',num2cell(handles.prmmap{handles.i,1})');
    set(handles.table_thresh,'Data',num2cell(handles.thresh));
    guidata(gcbo,handles);
end

% --- Executes on button press in button_crop_add.
function button_crop_add_Callback(~, ~, handles)
if isempty(handles.cutoff)
    handles.cutoff = [1,0,1];
else
    handles.cutoff(end+1,:) = handles.cutoff(end,:);
end
set(handles.table_crop,'Data',num2cell(handles.cutoff));
guidata(gcbo,handles);

% --- Executes on button press in button_crop_del.
function button_crop_del_Callback(~, ~, handles)
if ~isempty(handles.cuti)
    handles.cutoff(end,:) = [];
    set(handles.table_crop,'Data',num2cell(handles.cutoff));
    guidata(gcbo,handles);
end

% --- Executes on button press in button_C_add.
function button_C_add_Callback(~, ~, handles)
handles.prmmap{handles.i,1}(end+1,:) = handles.prmmap{handles.i,1}(end,:);
set(handles.table_C,'Data',num2cell(handles.prmmap{handles.i,1})');
guidata(gcbo,handles);

% --- Executes on button press in button_C_del.
function button_C_del_Callback(~, ~, handles)
if ~isempty(handles.Ci)
    handles.prmmap{handles.i,1}(handles.Ci,:) = [];
    set(handles.table_C,'Data',num2cell(handles.prmmap{handles.i,1})');
    guidata(gcbo,handles);
end

% --- Saves PRM options to a new default file.
function button_save_Callback(~, ~, handles)
if ~isempty(handles.prmmap)
    fpath = fullfile(fileparts(which('cmi')),'PRMdefs');
    dval = get(handles.popup_defs,'Value');
    if dval>1
        str = handles.defs(dval-1).name;
    else
        str = '';
    end
        
    str = inputdlg('PRM model name:','Save PRM Model',1,{str});
    str{1}(str{1}==' ') = []; % remove spaces
    if ~isempty(str{1})
        fpath = fullfile(fpath,['PRMdef_',str{1},'.mat']);
        answer = 'Yes';
        if exist(fpath,'file')
            answer = questdlg('File already exists! Replace?');
        end
        if strcmp(answer,'Yes')
            name = str{1};
            thresh = handles.thresh;
            cutoff = handles.cutoff;
            cmap = handles.cmap;
            prmmap = handles.prmmap;
            SPopts = handles.SPopts;
            save(fpath,'name','thresh','cutoff','cmap','prmmap','SPopts');
        end
    end
end

function edit_maxnp_Callback(~, ~, handles)
n = str2double(get(handles.edit_maxnp,'String'));
if ~isnan(n) && (n>0)
    handles.SPopts.Nmax = round(n);
end
set(handles.edit_maxnp,'String',num2str(handles.SPopts.Nmax));
guidata(gcbo,handles);

% --- Executes on button press in checkbox_showScatter.
function checkbox_showScatter_Callback(hObject, ~, handles)
handles.SPopts.show = logical(get(hObject,'Value'));
guidata(gcbo,handles);


% --- Executes on button press in button_moveUp.
function button_moveUp_Callback(~, ~, handles)
sel = get(handles.list_regions,'Value');
if sel>1
    sel = [sel-1,sel];
    handles.prmmap(sel,:) = handles.prmmap(sel([2,1]),:);
    handles.cmap(sel,:) = handles.cmap(sel([2,1]),:);
    set(handles.list_regions,'Value',sel(1),'String',handles.prmmap(:,2));
    guidata(gcbo,handles);
end

% --- Executes on button press in button_moveDwn.
function button_moveDwn_Callback(~, ~, handles)
sel = get(handles.list_regions,'Value');
if sel<size(handles.prmmap,1)
    sel = [sel,sel+1];
    handles.prmmap(sel,:) = handles.prmmap(sel([2,1]),:);
    handles.cmap(sel,:) = handles.cmap(sel([2,1]),:);
    set(handles.list_regions,'Value',sel(2),'String',handles.prmmap(:,2));
    guidata(gcbo,handles);
end

% Set scatterplot X-vector
function edit_SPdimX_Callback(hObject, ~, handles)
val = str2double(get(hObject,'String'));
if isnan(val)
    set(hObject,'String','');
end
handles.SPopts.Xvec = val;
guidata(gcbo,handles);

% Set scatterplot Y-vector
function edit_SPdimY_Callback(hObject, ~, handles)
val = str2double(get(hObject,'String'));
if isnan(val)
    set(hObject,'String','');
end
handles.SPopts.Yvec = val;
guidata(gcbo,handles);

% Set scatterplot X-min
function edit_SPminX_Callback(hObject, ~, handles)
val = str2double(get(hObject,'String'));
if isnan(val) || isnan(handles.SPopts.Xmin) || (val<handles.SPopts.Xmin)
    if isnan(val)
        set(hObject,'String','');
    end
    handles.SPopts.Xmin = val;
    guidata(gcbo,handles);
else
    if isnan(handles.SPopts.Xmin)
        str = '';
    else
        str = num2str(handles.SPopts.Xmin);
    end
    set(hObject,'String',str);
end

% Set scatterplot Y-min
function edit_SPminY_Callback(hObject, ~, handles)
val = str2double(get(hObject,'String'));
if isnan(val) || isnan(handles.SPopts.Ymin) || (val<handles.SPopts.Ymin)
    if isnan(val)
        set(hObject,'String','');
    end
    handles.SPopts.Ymin = val;
    guidata(gcbo,handles);
else
    if isnan(handles.SPopts.Ymin)
        str = '';
    else
        str = num2str(handles.SPopts.Ymin);
    end
    set(hObject,'String',str);
end

% Set scatterplot Y-max
function edit_SPmaxY_Callback(hObject, ~, handles)
val = str2double(get(hObject,'String'));
if isnan(val) || isnan(handles.SPopts.Ymin) || (val>handles.SPopts.Ymin)
    if isnan(val)
        set(hObject,'String','');
    end
    handles.SPopts.Ymax = val;
    guidata(gcbo,handles);
else
    if isnan(handles.SPopts.Ymax)
        str = '';
    else
        str = num2str(handles.SPopts.Ymax);
    end
    set(hObject,'String',str);
end

% Set scatterplot X-max
function edit_SPmaxX_Callback(hObject, ~, handles)
val = str2double(get(hObject,'String'));
if isnan(val) || isnan(handles.SPopts.Xmin) || (val>handles.SPopts.Xmin)
    if isnan(val)
        set(hObject,'String','');
    end
    handles.SPopts.Xmax = val;
    guidata(gcbo,handles);
else
    if isnan(handles.SPopts.Xmax)
        str = '';
    else
        str = num2str(handles.SPopts.Xmax);
    end
    set(hObject,'String',str);
end

% --- Executes on button press in checkbox_filtChk.
function checkbox_filtChk_Callback(hObject, ~, handles)
handles.filtchk = logical(get(hObject,'Value'));
guidata(gcbo,handles);

% Select pre-filter type
function popup_filt_Callback(hObject,~,handles)
str = get(hObject,'String');
handles.filttype = str{get(hObject,'Value')};
guidata(gcbo,handles);

% Size of filter window/strength
function edit_filtSize_Callback(hObject, ~, handles)
val = sscanf(get(hObject,'String'),'%f %f %f');
nv = length(val);
if any(nv==(1:2))
    if nv==1
        val = [val,val];
    end
    handles.filtstr = val;
    guidata(gcbo,handles);
else
    warning('Invalid filter strength.')
    set(hObject,'String',sprintf('% .1f % .1f',handles.filtstr));
end
