% RegClass function
function setTchk(self,hObject,~)
% Select whether or not to apply initial transform in Elastix call

val = [];
if isa(hObject,'matlab.ui.control.UIControl') && strcmp(hObject.Tag,'checkbox_useExistingT')
    val = hObject.Value;
elseif islogical(hObject)
    val = hObject;
end
if ~isempty(val)
    self.elxObj.setTx0par('check',logical(val));
    set(self.h.checkbox_useExistingT,'Value',val,'Enable','on');
end
