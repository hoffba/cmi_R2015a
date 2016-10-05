function varargout = multiTformGUI(varargin)
% MULTITFORMGUI MATLAB code for multiTformGUI.fig
%      MULTITFORMGUI, by itself, creates a new MULTITFORMGUI or raises the existing
%      singleton*.
%
%      H = MULTITFORMGUI returns the handle to a new MULTITFORMGUI or the handle to
%      the existing singleton*.
%
%      MULTITFORMGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MULTITFORMGUI.M with the given input arguments.
%
%      MULTITFORMGUI('Property','Value',...) creates a new MULTITFORMGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before multiTformGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to multiTformGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @multiTformGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @multiTformGUI_OutputFcn, ...
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

% --- Executes just before multiTformGUI is made visible.
function multiTformGUI_OpeningFcn(hObject, ~, handles, varargin)

% Choose default command line output for multiTformGUI
handles.output = hObject;
handles.table_TF.Data = {};
handles.table_im.Data = {};
setappdata(hObject,'tp',struct('chain',{},'fname',{},'im',{},'jac',{}));
setappdata(hObject,'tsel',[]);
setappdata(hObject,'isel',[]);
setappdata(hObject,'odir','');

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = multiTformGUI_OutputFcn(~, ~, handles) 
varargout{1} = handles.output;

function button_TFadd_Callback(~, ~, handles)
[fname,fdir] = uigetfile('*.txt','TransformParameters');
if ischar(fname)
    tp = getappdata(handles.figure1,'tp');
    i = length(tp)+1;
    tp(i).fname = fullfile(fdir,fname);
    tp(i).im = {};
    tp(i).jac = false;
    if i==1
        tp(i).chain = 1;
        tval = true;
    else
        tp(i).chain = tp(i-1).chain;
        tval = false;
    end
    
    % Update table
    dind = max(1,length(tp(i).fname)-100);
    handles.table_TF.Data = [handles.table_TF.Data;
        {tval,['... ',tp(i).fname(dind:end)],tp(i).jac}];
    
    setappdata(handles.figure1,'tp',tp);
end

function button_TFremove_Callback(~, ~, handles)
i = getappdata(handles.figure1,'tsel');
if i
    tp = getappdata(handles.figure1,'tp');
    tp(i) = [];
    setappdata(handles.figure1,'tp',tp);
    handles.table_TF.Data(i,:) = [];
end

function button_Iadd_Callback(~, ~, handles)
i = getappdata(handles.figure1,'tsel');
if i
    [fname,fdir] = uigetfile('*.mhd','Images','MultiSelect','on');
    if ischar(fdir)
        if ischar(fname)
            fname = {fname};
        end
        % Determine if interp should be nearest neighbor:
        nn = cellfun(@(x)~isempty(regexp(x,'VOI|ADC','Once')),fname)';
        tp = getappdata(handles.figure1,'tp');
        %j = size(tp(i).im,1)+(1:length(fname));
        % Set output file names:
        oname = cellfun(@(x)[x(1:end-4),'_mR',x(end-3:end)],fname,'UniformOutput',false)';
        im = [tp(i).im;num2cell(nn),...
              fullfile(fdir,fname)',...cellfun(@(x)fullfile(fdir,x),fname,'UniformOutput',false)',...
              oname];%cellfun(@(x){sprintf('TransformedImage.%02u.%02u.%02u.mhd',tp(i).chain,i,x)},num2cell(j))'];
        tp(i).im = im;
        setappdata(handles.figure1,'tp',tp);
        if ~isempty(im)
            im(:,2) = cellfun(@(x)['... ',x(max(1,length(x)-80):end)],im(:,2),'UniformOutput',false);
        end
        set(handles.table_im,'Data',im);
    end
end

function button_Iremove_Callback(~, ~, handles)
i = getappdata(handles.figure1,'tsel');
j = getappdata(handles.figure1,'isel');
if i
    tp = getappdata(handles.figure1,'tp');
    tp(i).im(j,:) = [];
    setappdata(handles.figure1,'tp',tp);
    im = tp(i).im;
    if ~isempty(im)
        im(:,2) = cellfun(@(x)['... ',x(max(1,length(x)-80):end)],im(:,2),'UniformOutput',false);
    end
    set(handles.table_im,'Data',im);
end

function button_TFup_Callback(~, ~, handles)
i = getappdata(handles.figure1,'tsel');
if i
    tp = getappdata(handles.figure1,'tp');
    tp([i-1,i]) = tp([i,i-1]);
    setappdata(handles.figure1,'tp',tp);
    handles.table_TF.Data([i-1,i],:) = handles.table_TF.Data([i,i-1],:);
    pause(0.1);
    setappdata(handles.figure1,'tsel',i-1);
end

function button_TFdown_Callback(~, ~, handles)
i = getappdata(handles.figure1,'tsel');
if i
    tp = getappdata(handles.figure1,'tp');
    tp([i,i+1]) = tp([i+1,i]);
    setappdata(handles.figure1,'tp',tp);
    handles.table_TF.Data([i,i+1],:) = handles.table_TF.Data([i+1,i],:);
    pause(0.1);
    setappdata(handles.figure1,'tsel',i+1);
end

function button_clear_Callback(~, ~, handles)
setappdata(handles.figure1,'tp',struct('fin',{},'fname',{},'im',{}));
setappdata(handles.figure1,'tsel',0);
setappdata(handles.figure1,'isel',0);
handles.table_TF.Data = {};
handles.table_im.Data = {};

function button_help_Callback(~, ~, ~)
str = {'TRANSFORMS:';...
       ['Add TransformParameters.*.txt files in forward order ',...
        '(transforms are applied backwards, from end of list).',...
        'Mark final transformations in the chain in the first column. ',...
        'Images at the end will be transformed up through the chain until ',...
        'a final transform is reached.'];...
       '';...
       'IMAGES:';...
       'MHD images may be added to any step in the chain to be transformed from that point upward.';...
       '';...
       'OUTPUT:';...
       ['All transforms, images, and a .log file are saved in the same selected directory. , '...
        'Images have the following format: [fname].[chain#].[step#].[image#].mhd ',...
        'and original file names are saved in the log.']};
msgbox(str,'MultiTransform HELP','help','modal');

function button_odir_Callback(~, ~, handles)
odir = getappdata(handles.figure1,'odir');
if isempty(odir)
    odir = pwd;
else
    odir = fileparts(odir);
end
odir = uigetdir(odir,'Output directory');
if odir~=0
    setappdata(handles.figure1,'odir',odir);
    set(handles.edit_odir,'String',odir);
end

function edit_odir_Callback(hObject, ~, handles)
odir = hObject.String;
stat = false;
if ~exist(odir,'dir')
    answer = questdlg('Directory not found. Create new?','MkDir','Yes','Cancel','Cancel');
    if strcmp(answer,'Yes')
        stat = mkdir(odir);
    else
        set(handles.edit_odir,'String',getappdata(handles.figure1,'odir'));
    end
end
if stat
    setappdata(handles.figure1,'odir',odir);
end

% --- Executes when selected cell(s) is changed in table_TF.
function table_TF_CellSelectionCallback(~, eventdata, handles)
i = eventdata.Indices;
if ~isempty(i)
    i = i(1);
    tp = getappdata(handles.figure1,'tp');
    im = tp(i).im;
    if ~isempty(im)
        im(:,2) = cellfun(@(x)['... ',x(max(1,length(x)-80):end)],im(:,2),'UniformOutput',false);
    end
    set(handles.table_im,'Data',im);
end
setappdata(handles.figure1,'tsel',i);

function table_TF_CellEditCallback(hObject, eventdata, handles)
row = eventdata.Indices(1);
col = eventdata.Indices(2);
tp = getappdata(handles.figure1,'tp');
if col==1
    if row==1
        % First transform MUST be a final
        handles.table_TF.Data{1} = true;
    else
        chain = cumsum([true,hObject.Data{2:end,1}]);
        for row = 1:length(tp)
            tp(row).chain = chain(row);
        end
    end
elseif col==3
    tp(row).jac = eventdata.NewData;
end
setappdata(handles.figure1,'tp',tp);

function table_im_CellEditCallback(hObject, eventdata, handles)
tp = getappdata(handles.figure1,'tp');
i = getappdata(handles.figure1,'tsel');
ind = eventdata.Indices;
if ind(2)==1 % NN check
    tp(i).im{ind(1),1} = eventdata.NewData;
elseif ind(2)==3 % Output name
    str = eventdata.NewData;
    if ~strcmp(str(end-3:end),'.mhd')
        str = [str,'.mhd'];
    end
    tp(i).im{ind(1),3} = str;
    hObject.Data{ind(1),3} = str;
end
setappdata(handles.figure1,'tp',tp);

function table_im_CellSelectionCallback(~, eventdata, handles)
i = eventdata.Indices;
if ~isempty(i)
    i = i(1);
end
setappdata(handles.figure1,'isel',i);

function button_START_Callback(~, ~, handles)
multiTform(getappdata(handles.figure1,'tp'),...
    getappdata(handles.figure1,'odir'));

    
