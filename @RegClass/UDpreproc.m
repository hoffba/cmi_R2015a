% RegClass function
function UDpreproc(self,x,val)
% Handle Callbacks from 

if ishandle(x)
    tag = get(x,'Tag');
    switch tag
        case 'popup_filterType'
            str = get(x,'String');
            self.ftype = str{get(x,'Value')};
        case {'edit_preFilter','edit_VOIdilate'}
            str = {'filtN','dilateN'};
            fldni = strcmp(tag,'edit_VOIdilate')+1;
            i = str2num(get(x,'String'));
            if isempty(i) || (length(i)~=3) || any(isnan(i)|isinf(i))
                set(x,'String',num2str(self.(str{fldni})));
            else
                self.(str{fldni}) = i;
            end
        case {'edit_clampRef','edit_clampHom'}
            row = strcmp(tag,'edit_clampHom')+1;
            i = str2num(get(x,'String'));
            if isempty(i) || (length(i)~=2) || any(isnan(i))
                set(x,'String',num2str(self.clamp(row,:)));
            else
                self.clamp(row,:) = sort(i);
            end
        case 'checkbox_histeq'
            self.histeq = get(x,'Value');
        otherwise
            warning(['Unknown tag:',tag])
    end
elseif ischar(x) && isprop(self,x) && (nargin==3)
    switch x
        case 'ftype'
            h = self.h.popup_filterType;
            str = get(h,'String');
            i = find(strcmpi(val,str),1);
            if ischar(val) && ~isempty(i)
                set(h,'Value',i);
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
                if strcmp(x,'clampRef');
                    i = 1;
                    str = 'edit_clampRef';
                else
                    i = 2;
                    str = 'edit_clampHom';
                end
                set(self.h.(str),'String',num2str(val));
                self.clamp(i,:) = val;
            end
        case 'histeq'
            if islogical(val)
                set(self.h.checkbox_histeq,'Value',val);
                self.histeq = val;
            end
    end      
else
    warning('RegClass.UDpreproc() : Invalid inputs');
end