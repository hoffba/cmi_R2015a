% CMIclass function
function setDrawMode(self,x,~)
% Set draw mode function
%   0 = roipoly
%   1 = imfreehand

dm = self.drawMode;
if nargin==0
    dm = strcmp(questdlg('Draw Mode:','Draw Mode:',...
                         'roipoly','imfreehand','roipoly'),'imfreehand');
elseif (nargin==2) && islogical(x)
    dm = x;
elseif (nargin==3) && strcmp(x.Tag,'tools_voi_drawmode')
    dm = ~self.drawMode;
    if dm
        dstr = 'free';
    else
        dstr = 'poly';
    end
    x.Label = ['Draw: ',dstr];
end
self.drawMode = dm;