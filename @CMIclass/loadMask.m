% CMIclass function
% Load mask
function stat = loadMask(self,x,optsel)
if (nargin==1) || isa(x,'matlab.ui.container.Menu') || ~ischar(x)
    x = '';
    optsel = '';
end
if nargin<3
    optsel = '';
end
stat = self.img.loadMask(x,optsel);
disp(['VOI Volume : ',num2str(numel(self.img.mask.mat(self.img.mask.mat>0))*prod(self.img.voxsz)/1e6),' L']);
if stat && self.guicheck
    self.dispUDroi;
    self.dispUDhist;
end