% CMIclass function
function setVOImode(self,x,~)
% roifun specifies the add/replace/remove function applied to the current mask
%       = 'Auto' -> determines whether to add/remove based on what's selected
%       = 'Replace'
%       = 'Add'
%       = 'Intersect'
%       = 'Remove'
opts = {'Auto','Replace','Add','Intersect','Remove'};
str = [];
if (nargin>1)
    if (nargin==3) && ishandle(x) && strcmp(get(x,'Tag'),'popup_voimode')
        str = get(x,'String');
        str = str{get(x,'Value')};
    elseif ischar(x) && any(strcmp(x,opts))
        str = x;
    end
end
if ~isempty(str)
    self.roifun = str;
end