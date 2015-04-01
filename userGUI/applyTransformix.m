function varargout = applyTransformix(varargin)
% APPLYTRANSFORMIX MATLAB code for applyTransformix.fig
%      APPLYTRANSFORMIX, by itself, creates a new APPLYTRANSFORMIX or raises the existing
%      singleton*.
%
%      H = APPLYTRANSFORMIX returns the handle to a new APPLYTRANSFORMIX or the handle to
%      the existing singleton*.
%
%      APPLYTRANSFORMIX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in APPLYTRANSFORMIX.M with the given input arguments.
%
%      APPLYTRANSFORMIX('Property','Value',...) creates a new APPLYTRANSFORMIX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before applyTransformix_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to applyTransformix_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help applyTransformix

% Last Modified by GUIDE v2.5 04-Mar-2015 15:38:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @applyTransformix_OpeningFcn, ...
                   'gui_OutputFcn',  @applyTransformix_OutputFcn, ...
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

% --- Executes just before applyTransformix is made visible.
function applyTransformix_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to applyTransformix (see VARARGIN)

% Choose default command line output for applyTransformix
handles.output = hObject;
set(handles.uitable1,'Data',{});

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes applyTransformix wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = applyTransformix_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in button_tp.
function button_tp_Callback(~, ~, handles)
[fname,fpath] = uigetfile('*.txt','Select Transform File:');
if ischar(fname)
    set(handles.text_tp,'String',fullfile(fpath,fname));
end

% --- Executes on button press in button_odir.
function button_odir_Callback(~, ~, handles)
fpath = uigetdir(get(handles.text_odir,'String'),'Select Output Directory');
if iscellstr(fpath) || ischar(fpath)
    set(handles.text_odir,'String',fpath);
end

% --- Executes on button press in button_add.
function button_add_Callback(~, ~, handles)
[fname,fpath] = uigetfile('*.mhd','Select Images to Transform',...
                          'MultiSelect','on');
if iscellstr(fname) || ischar(fname)
    if ischar(fname)
        fname = {fname};
    end
    vchk = cellfun(@(x)~isempty(strfind(x,'VOI')),fname,'UniformOutput',false)';
    fpath = {fpath};
    nf = length(fname);
    if nf>1
        fpath = repmat(fpath,[nf,1]);
    end
    set(handles.uitable1,'Data',[get(handles.uitable1,'Data');...
                                 vchk,fname',fpath]);
end

% --- Executes on button press in button_remove.
function button_remove_Callback(~, ~, handles)
ind = inputdlg('What file #s to remove:','Remove Files');
ind = round(str2num(ind{1}));
if ~isempty(ind)
    tdata = get(handles.uitable1,'Data');
    ind((ind<1)|(ind>size(tdata,1))) = [];
    tdata(ind,:) = [];
    set(handles.uitable1,'Data',tdata);
end

% --- Executes on button press in button_start.
function button_start_Callback(~, ~, handles)
% Grab input options:
tdata = get(handles.uitable1,'Data');
nf = size(tdata,1);
chk = logical(cell2mat(get([handles.checkbox_jac,...
                            handles.checkbox_jacmat,...
                            handles.checkbox_def],'Value')));
tp = get(handles.text_tp,'String');
odir = get(handles.text_odir,'String');
if ~exist(tp,'file')
    errordlg('TransformParameter file does not exist.');
elseif ~isdir(odir)
    errordlg('Output directory does not exist.');
elseif (any(chk) || (nf>0))
    elxObj = ElxClass;
    vchk = logical(cell2mat(tdata(:,1)));
    if any(vchk)
        % Need to save copy of TransformParameter file with NearestNeighbor
        elxObj.loadTx(tp);
        elxObj.setTx0par('ResampleInterpolator','FinalNearestNeighborInterpolator');
        tpNN = [tp(1:end-4),'_NN.txt'];
        elxObj.saveTx0(tpNN);
    end
    ni = max(1,nf);
    str = cell(ni,1);
    for i = 1:max(1,nf)
        if vchk(i)
            tpfn = tpNN;
        else
            tpfn = tp;
        end
        if nf>0
            tfn = fullfile(tdata{i,3},tdata{i,2});
            inC = {'in',tfn,'outfn',[tfn(1:end-4),'_R.mhd']};
        end
        str(i) = elxObj.tfxCmd(odir,'tp',tpfn,'jac',chk(1),...
                    'jacmat',chk(2),'def',chk(3),inC{:});
        chk(:) = false;
    end
    % Start Transformix in new xterm window:
    namestr = '(Transformix) Applying Transform to Selected Images';
    system(['xterm -geometry 170x50 -T "',namestr,'"',...
            ' -e ''',strcat(str{:}),';csh''&']);
end
