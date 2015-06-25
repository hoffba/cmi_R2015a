% RegClass function
function UDtfx(self,varargin)
% GUI callback to update Transformix options
%   OR user-input Name/Value pairs:
%       'par'       parameter file name
%       'out'       output directory
%       'jac'       true/false - output 
%       'jacmat'
%       'def'
%       'nn'
%       'fnames'

C = varargin;
if (nargin==1)
    return;
end
if isa(varargin{1},'matlab.ui.control.UIControl')
    % Source of callback determines what is updated
    hObject = varargin{1};
    switch hObject.Tag
        case 'button_TfxPar'
            [fname,fpath] = uigetfile('*.txt','Select Transformix Parameters');
            if (fname~=0)
                % * Default Output to same directory
                fname = fullfile(fpath,fname);
                C = {'par',fname,'out',fpath};
                self.h.text_TfxPar.String = fname;
                self.h.text_TfxOut.String = fpath;
            end
        case 'button_TfxOut'
            fpath = uigetdir(pwd,'Select Output Directory');
            if (fpath~=0)
                C = {'out',fpath};
                self.h.text_TfxOut.String = fpath;
            end
        case 'checkbox_TfxJac'
            C = {'jac',logical(hObject.Value)};
        case 'checkbox_TfxJacmat'
            C = {'jacmat',logical(hObject.Value)};
        case 'checkbox_TfxDef'
            C = {'def',logical(hObject.Value)};
        case 'button_TfxAdd'
            [fnames,fpath] = uigetfile('*.mhd','Select Images to Transform',...
                'MultiSelect','on');
            if (fpath~=0)
                if ischar(fnames)
                    fnames = {fnames};
                end
                fnames = fnames(:);
                vchk = ~cellfun(@isempty,strfind(fnames,'VOI'));
                self.h.table_Tfx.Data = [self.h.table_Tfx.Data;...
                    num2cell(vchk),fnames,repmat({fpath},length(fnames),1)];
                vchk = [self.Tfx.nn ; vchk];
                fnames = [self.Tfx.fnames;...
                    cellfun(@(x)fullfile(fpath,x),fnames,'UniformOutput',false)];
                C = {'fnames',fnames,'nn',vchk};
            end
        case 'button_TfxRemove'
            if hObject.Value
                fcn = @self.UDtfx;
                hObject.ForegroundColor = 'red';
                hObject.FontWeight = 'bold';
                str = 'on';
            else
                fcn = '';
                hObject.ForegroundColor = 'black';
                hObject.FontWeight = 'Normal';
                str = 'off';
            end
            self.h.table_Tfx.CellSelectionCallback = fcn;
            self.h.table_Tfx.SelectionHighlight = str;
            C = {};
        case 'table_Tfx'
    end
elseif isa(varargin{1},'matlab.ui.control.Table')
    hObject = varargin{1};
    evnt = varargin{2};
    switch evnt.EventName
        case 'CellEdit'
            % Only the first column is editable (VOI checkboxes)
            C = {'nn',[hObject.Data{:,1}]};
        case 'CellSelection'
            % Only used to delete rows from table
            if isempty(evnt.Indices)
                C = {};
            else
                ind = evnt.Indices(1);
                tC = hObject.Data;
                tC(ind,:) = [];
                hObject.Data = tC;
                tnn = [tC{:,1}];
                if isempty(tnn)
                    tnn = false(0,1);
                end
                C = {'fnames',tC(:,2),'nn',tnn};
            end
    end
end
if ~isempty(C) && iscellstr(C(1:2:end))
    p = inputParser;
    p.CaseSensitive = false;
    addParameter(p,'par',self.Tfx.par,@(x)ischar(x)&&exist(x,'file'));
    addParameter(p,'out',self.Tfx.out,@(x)ischar(x)&&isdir(x));
    addParameter(p,'jac',self.Tfx.jac,@(x)islogical(x));
    addParameter(p,'jacmat',self.Tfx.jacmat,@(x)islogical(x));
    addParameter(p,'def',self.Tfx.def,@(x)islogical(x));
    addParameter(p,'nn',self.Tfx.nn,@(x)islogical(x));
    addParameter(p,'fnames',self.Tfx.fnames,@(x)iscellstr(x));
    parse(p,C{:});
    instr = setxor(p.Parameters,p.UsingDefaults);
    
    % Make sure fnames and nn have equal length
    if any(ismember(instr,{'fnames','nn'}))
        nf = length(p.Results.fnames);
        nnn = length(p.Results.nn);
        if nf>nnn
            p.Results.nn(nf+1:end) = false;
        elseif nnn>nf
            p.Results.nn(nf+1:end) = [];
        end
    end
    for i = 1:length(instr)
        self.Tfx.(instr{i}) = p.Results.(instr{i});
    end
end

