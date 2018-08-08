% CMIclass function
% Set colormap for selected image (background or overlay)
function setCMap(self,x,~)
opts = {'gray','jet','hsv','hot','cool','spring','summer','autumn','winter',...
        'bone','copper','pink','lines'};
str = [];
if (nargin>1)
    if (nargin==3) && ishandle(x) && strcmp(get(x,'Tag'),'popup_colormap')
        str = get(x,'String');
        str = str{get(x,'Value')};
    elseif ischar(x) && any(strcmp(x,opts))
        str = x;
    end
end
if ~isempty(str)
    if self.guicheck
        set(self.h.popup_colormap,'Value',find(strcmp(str,opts),1));
    end
    if self.bgcheck
        self.bgcmap = str;
    else
        self.cmap = str;
    end
    self.dispUDcmap;
end
