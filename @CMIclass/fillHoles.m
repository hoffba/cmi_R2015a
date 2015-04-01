% CMIclass function
% Fills holes in current mask
function fillHoles(self,~,~)
self.img.mask.fillHoles;
self.dispUDroi;