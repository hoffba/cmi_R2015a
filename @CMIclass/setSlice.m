% CMIclass function
% Set slice in current view
function setSlice(self,x,~)
val = [];
% Determine input value:
if nargin>1
    if (nargin==3) && ishghandle(x)
        % input from GUI
        if strcmp(x.Tag,'slider_slice')
                val = x.Value;
        elseif strcmp(x.Tag,'edit_slc')
                val = str2double(x.String);
        end
    elseif isnumeric(x)
        val = x;
    end
end
val = round(val);
% If valid, set properties
if ((val > 0) && (val <= self.img.dims(self.orient)))
    % Update GUI properties
    if self.guicheck
        self.h.slider_slice.Value = val;
        self.h.edit_slc.String = num2str(val);
    end
    % Update CMIclass properties and update display
    self.slc(self.orient) = val;
    self.dispUDslice;
    self.dispUDhist(1);
end