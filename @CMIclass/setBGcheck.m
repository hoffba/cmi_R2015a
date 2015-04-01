% CMIclass function
% Set background check
function setBGcheck(self,x,~)
val = [];
if (nargin>1)
    if (nargin==3) && ishandle(x) && strcmp(get(x,'Tag'),'bgcheck')
        val = logical(get(x,'Value'));
    elseif isnumeric(x) && ~isempty(x)
        val = logica(x(1));
    end
end
if ~isempty(val)
    self.bgcheck = val;
    % Update GUI colormap options
    if isstruct(self.h)
        if val
            tcmap = self.bgcmap;
        else
            tcmap = self.cmap;
        end
        [tmin,tmax] = self.getColorLimits;
        [imin,imax] = self.getColorMinMax;
        cpad = (imax - imin) / 2;
        set(self.h.slider_cmin,'Min',(imin-cpad),'Max',(imax+cpad),...
                                'Value',tmin,'Enable','On');
        set(self.h.slider_cmax,'Min',(imin-cpad),'Max',(imax+cpad),...
                                'Value',tmax,'Enable','On');
        set(self.h.edit_cmin,'String',num2str(tmin),'Enable','On');
        set(self.h.edit_cmax,'String',num2str(tmax),'Enable','On');
        tval = find(strcmp(tcmap,get(self.h.popup_colormap,'String')),1);
        set(self.h.popup_colormap,'Value',tval);
    end
end