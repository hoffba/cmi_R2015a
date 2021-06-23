% CMIclass function
% Set image 4D label
function setLabel(self,x,~)
str = [];
if (nargin==3) && ishandle(x) && strcmp(get(x,'Tag'),'edit_veclabel')
    str = get(x,'String');
elseif self.img.check && ischar(x)
    str = x;
end
if ~isempty(str)
    self.img.setLabel(self.vec,str);
end