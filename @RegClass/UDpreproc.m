% RegClass function
function UDpreproc(self,varargin)
% Handle Callbacks from 

if self.guicheck && (nargin>1) && ishandle(varargin{1}(1))
    h = varargin{1}(1);
    tag = get(h,'Tag');
    switch tag
        case 'popup_filterType'
            str = get(h,'String');
            self.ftype = str{get(h,'Value')};
        case 'checkbox_histeq'
            self.histeq = get(h,'Value');
        case 'table_opts'
            irow = varargin{2}.Indices(1);
            icol = varargin{2}.Indices(2);
            val = varargin{2}.NewData;
            oldval = varargin{2}.PreviousData;
            switch irow
                case 1 % min
                    if isscalar(val) && (val<self.clamp(icol,2))
                        self.clamp(icol,1) = val;
                    else
                        h.Data{irow,icol} = oldval;
                    end
                case 2 % max
                    if isscalar(val) && (val>self.clamp(icol,1))
                        self.clamp(icol,2) = val;
                    else
                        h.Data{irow,icol} = oldval;
                    end
                case {3,4,5} % filtX, filtY, filtZ
                    if isscalar(val) && (val>=0)
                        self.filtN(icol,irow-2) = val;
                    else
                        h.Data{irow,icol} = oldval;
                    end
                case {6,7,8} % DilateX, DilateY, DilateZ
                    if isscalar(val) && (val>=0)
                        self.dilateN(icol,irow-5) = val;
                    else
                        h.Data{irow,icol} = oldval;
                    end
                case 9 % UnMaskVal
                    if isscalar(val) && ~isinf(val)
                        self.unmaskval(icol) = val;
                    else
                        h.Data{irow,icol} = oldval;
                    end
            end
        otherwise
            warning(['Unknown tag:',tag])
    end
elseif (nargin>=3) && ischar(varargin{1})
    % Parse name/value pairs
    f = varargin(1:2:end);
    v = varargin(2:2:end);
    if length(f)~=length(v)
        error('RegClass/UDpreproc -> Invalid input. Requires Name/Value pairs.')
    end
    for i = 1:length(f)
        val = v{i};
        switch lower(f{i})
            case 'ftype'
                h = self.h.popup_filterType;
                str = get(h,'String');
                ind = find(strcmpi(val,str),1);
                if ischar(val) && ~isempty(i)
                    set(h,'Value',ind);
                    self.ftype = val;
                end
            case 'filtn'
                if isnumeric(val) && (size(val,1)==2) && (size(val,2)==3) && ~any(isnan(val(:))|isinf(val(:)))
                    if self.guicheck
                        self.h.table_opts.Data(3:5,:) = num2cell(val');
                    end
                    self.filtN = val;
                end
            case 'dilaten'
                if isnumeric(val) && (size(val,1)==2) && (size(val,2)==3) && ~any(isnan(val(:))|isinf(val(:)))
                    if self.guicheck
                        self.h.table_opts.Data(6:8,:) = num2cell(val');
                    end
                    self.dilateN = val;
                end
            case 'clamp'
                if isnumeric(val) && (size(val,1)==2) && (size(val,2)==2) && ~any(isnan(val(:)))
                    val = sort(val,2);
                    if self.guicheck
                        self.h.table_opts.Data(1:2,:) = num2cell(val');
                    end
                    self.clamp = val;
                end
            case 'histeq'
                if islogical(val)
                    if self.guicheck
                        set(self.h.checkbox_histeq,'Value',val);
                    end
                    self.histeq = val;
                end
            case 'unmaskval'
                if isnumeric(val) && (numel(val)==2) && ~any(isinf(val))
                    if self.guicheck
                        self.h.table_opts.Data(9,:) = num2cell(val);
                    end
                    self.unmaskval = val;
                end
            otherwise
                warning('Property "%s" not set',f{i});
        end      
    end
else
    warning('RegClass.UDpreproc() : Invalid inputs');
end