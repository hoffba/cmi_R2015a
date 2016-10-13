% RegClass function
function UDpreproc(self,varargin)
% Handle Callbacks from 

if (nargin>1) && ishandle(varargin{1}(1))
    h = varargin{1}(1);
    tag = get(h,'Tag');
    switch tag
        case 'popup_filterType'
            str = get(h,'String');
            self.ftype = str{get(h,'Value')};
        case {'edit_preFilter','edit_VOIdilate'}
            str = {'filtN','dilateN'};
            fldni = strcmp(tag,'edit_VOIdilate')+1;
            i = str2num(get(h,'String'));
            if isempty(i) || (length(i)~=3) || any(isnan(i)|isinf(i))
                set(h,'String',num2str(self.(str{fldni})));
            else
                self.(str{fldni}) = i;
            end
        case {'edit_clampRef','edit_clampHom'}
            row = strcmp(tag,'edit_clampHom')+1;
            i = str2num(get(h,'String'));
            if isempty(i) || (length(i)~=2) || any(isnan(i))
                set(x,'String',num2str(self.clamp(row,:)));
            else
                self.clamp(row,:) = sort(i);
            end
        case 'checkbox_histeq'
            self.histeq = get(h,'Value');
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
        switch f{i}
            case 'ftype'
                h = self.h.popup_filterType;
                str = get(h,'String');
                ind = find(strcmpi(val,str),1);
                if ischar(val) && ~isempty(i)
                    set(h,'Value',ind);
                    self.ftype = val;
                end
            case 'filtN'
                if isnumeric(val) && (length(val)==3) && ~any(isnan(val)|isinf(val))
                    set(self.h.edit_preFilter,'String',num2str(val));
                    self.filtN = val;
                end
            case 'dilateN'
                if isnumeric(val) && (length(val)==3) && ~any(isnan(val)|isinf(val))
                    set(self.h.edit_VOIdilate,'String',num2str(val));
                    self.dilateN = val;
                end
            case {'clampRef','clampHom'}
                if isnumeric(val) && (length(val)==2)
                    val = sort(val); % Make sure it's in ascending order
                    if strcmp(f{i},'clampRef');
                        ind = 1;
                        str = 'edit_clampRef';
                    else
                        ind = 2;
                        str = 'edit_clampHom';
                    end
                    set(self.h.(str),'String',num2str(val));
                    self.clamp(ind,:) = val;
                end
            case 'histeq'
                if islogical(val)
                    set(self.h.checkbox_histeq,'Value',val);
                    self.histeq = val;
                end
            otherwise
                warning('Property "%s" not set',f{i});
        end      
    end
else
    warning('RegClass.UDpreproc() : Invalid inputs');
end