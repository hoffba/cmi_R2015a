% CMIclass function
% Set image values inside/outside mask to specific value
function setMaskVal(self,vec,mval,ival)

stat = false;
if (nargin<4)
    % Ask for use input:
    answer = inputdlg({'Images:','Mask Value (0/1):','New Image Value:'},...
        'SetMaskValue',1,{num2str(1:self.img.dims(4)),'0','0'});
    if ~isempty(answer)
        vec = str2num(answer{1});
        mval = str2double(answer{2});
        ival = str2double(answer{3});
        stat = true;
    end
else
    stat = true;
end
if stat
    stat = self.img.setMaskVal(vec,mval,ival);
    if stat
        self.dispUDslice;
    end
end