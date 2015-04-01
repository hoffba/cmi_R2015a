% CMIclass function
% Load mask
function stat = loadMask(self,x,~)
if (nargin==1) || (isnumeric(x) && ishandle(x)) || ~ischar(x)
    x = '';
end
stat = self.img.loadMask(x);
if stat
    self.dispUDroi;
    self.dispUDhist;
end