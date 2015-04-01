% CMIclass function
function connToggle(self,hObject,~)
% Toggles between 2D and 3D (default) connectivity space

if self.chk2D
    self.chk2D = false;
    str = '3D';
else
    self.chk2D = true;
    str = '2D';
end
if (nargin>1) && strcmp(get(hObject,'Tag'),'tools_connToggle')
    set(hObject,'Label',['Connectivity : ',str]);
end