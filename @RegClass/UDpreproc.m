% RegClass function
function UDpreproc(self,hObject,~)
% Handle Callbacks from 

if ishandle(hObject)
    tag = get(hObject,'Tag');
    switch tag
        case 'popup_filterType'
            str = get(hObject,'String');
            self.ftype = str{get(hObject,'Value')};
        case {'edit_preFilter','edit_VOIdilate'}
            str = {'filtN','dilateN'};
            fldni = strcmp(tag,'edit_VOIdilate')+1;
            i = str2num(get(hObject,'String'));
            if isempty(i) || (length(i)~=3) || any(isnan(i)|isinf(i))
                set(hObject,'String',num2str(self.(str{fldni})));
            else
                self.(str{fldni}) = i;
            end
        case {'edit_clampRef','edit_clampHom'}
            row = strcmp(tag,'edit_clampHom')+1;
            i = str2num(get(hObject,'String'));
            if isempty(i) || (length(i)~=2) || any(isnan(i))
                set(hObject,'String',num2str(self.clamp(row,:)));
            else
                self.clamp(row,:) = sort(i);
            end
        otherwise
            warning(['Unknown tag:',tag])
    end
else
    warning('RegClass.UDpreproc() : Invalid inputs');
end