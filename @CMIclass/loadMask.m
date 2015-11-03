% CMIclass function
% Load mask
function stat = loadMask(self,x,~)
if (nargin==1) || (isnumeric(x) && ishandle(x)) || ~ischar(x)
    x = '';
end
stat = self.img.loadMask(x);
['VOI Volume : ',num2str(numel(self.img.mask.mat(self.img.mask.mat>0))*prod(self.img.voxsz)/1e6),' L']
if stat
    self.dispUDroi;
    self.dispUDhist;
end