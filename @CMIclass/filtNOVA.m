% CMIclass function
function filtNOVA(self,~,~)
if self.img.check
    stat = self.img.filtNOVA;
    if stat
        self.dispUDimg;
    end
end