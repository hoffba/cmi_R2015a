% CMIclass function
function stat = saveMask(self,~,~)
if self.img.mask.check
    stat = self.img.saveMask;
end