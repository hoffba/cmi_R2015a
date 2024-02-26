% CMIclass function
% Set colormap for selected image (background or overlay)
function setCMap(self,x,~)
opts = {'jet','hsv','hot','cool','spring','summer','autumn','winter',...
        'gray','bone','copper','pink','lines','AtlasCTA'};
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
    if strcmp(str,'AtlasCTA')
        ls = linspace(0.25,1,25)';
        z = zeros(25,1);
        str = [ 0 ,  0 ,  0  ;...
                z ,  ls , z  ;... % green
                z ,  ls , ls ;... % cyan
                ls , z ,  z  ;... % red
                ls , ls , z  ;... % yellow
                ls , z ,  ls      % magenta
              ];
        self.ncolors = size(str,1);
        % Need to set min/max
        self.setClim([self.vec,91,599]);
    else
        self.ncolors = 128;
    end
    if self.bgcheck
        self.bgcmap = str;
    else
        self.cmap = str;
    end
    self.dispUDcmap;
end
