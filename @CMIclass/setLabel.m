% CMIclass function
% Set image 4D label
function setLabel(self,x,cstr)
if (nargin==3) && isa(x(1),'matlab.ui.control.UIControl') && strcmp(get(x,'Tag'),'edit_veclabel')
    cstr = {get(x,'String')};
    x = self.vec;
elseif self.img.check && ischar(x)
    cstr = {x};
    x = self.vec;
end
if ~isempty(cstr)
    self.img.setLabel(x,cstr);
end