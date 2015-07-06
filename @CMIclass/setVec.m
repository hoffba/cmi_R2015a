% CMIclass function
% Set 4D vector image index
function setVec(self,x,~)
val = [];
% Determine input value:
if (nargin>1)
    if (nargin==3) && ishghandle(x)
        if strcmp(x.Tag,'slider_4D')
            val = x.Value;
        elseif strcmp(x.Tag,'edit_vec')
            val = str2double(x.String);
        end
    elseif isnumeric(x)
        val = x;
    end
end
val = round(val);
% If valid, set properties
if ((val > 0) && (val <= self.img.dims(4)))
    % Update CMIclass properties and update display
    self.vec = val;
    if self.prmcheck % PRM needs re-calculation
        self.img.calcPRM(val);
    end
    % Update GUI properties
    if self.guicheck
        self.GUIupdate;
%         self.h.edit_veclabel.String = self.img.labels{val};
%         self.h.slider_4D.Value = val;
%         self.h.edit_vec.String = num2str(val);
%         [tmin,tmax] = self.getColorLimits;
%         [imin,imax] = self.getColorMinMax;
%         cpad = (imax - imin) / 2;
%         set(self.h.slider_cmin,'Min',(imin-cpad),'Max',(imax+cpad),...
%                                 'Value',tmin);
%         set(self.h.slider_cmax,'Min',(imin-cpad),'Max',(imax+cpad),...
%                                 'Value',tmax);
%         self.h.edit_cmin.String = num2str(tmin);
%         self.h.edit_cmax.String = num2str(tmax);
    end
%     self.dispUDslice;
    self.dispUDhist;
end
