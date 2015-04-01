% CMIclass function
% Filter image
function imgFilt(self,~,~)
if self.img.check
    self.img.imgFilt(self.vec);
    self.dispUDslice;
end