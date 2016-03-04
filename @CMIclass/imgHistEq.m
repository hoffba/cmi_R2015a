% CMIclass function
function imgHistEq(self,~,~)

stat = self.img.imgHistEq(self.vec);
if stat
    self.dispUDslice;
end