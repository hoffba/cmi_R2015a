% RegClass function
function setTchk(self,hObject,~)
% Select whether or not to apply initial transform in Elastix call

val = [];
if nargin==1
    val = ~self.elxObj.T0check; % toggle
elseif isa(hObject,'matlab.ui.control.UIControl') && strcmp(hObject.Tag,'checkbox_useExistingT')
    val = hObject.Value;
elseif islogical(hObject)
    val = hObject;
end
if ~isempty(val)
    self.elxObj.setTx0par('check',logical(val));
    if self.guicheck
        set(self.h.checkbox_useExistingT,'Value',val,'Enable','on');
    end
end
