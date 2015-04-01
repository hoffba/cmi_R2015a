% CMIcalss function
% Set PRM model
function setPRMmodel(self,~,~)
self.img.prm.setModel; % no inputs -> prompt user for input
if self.img.prm.check && self.prmcheck
    self.img.calcPRM(self.vec);
    self.dispUDcmap;
    self.dispUDimg;
end