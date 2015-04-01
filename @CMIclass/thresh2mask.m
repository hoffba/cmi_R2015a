% CMIclass function
% Add threshold to mask
function thresh2mask(self,~,~)
if self.img.check
    self.img.thresh2mask(self.vec);
    self.dispUDmask;
end