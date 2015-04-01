% CMIclass function
% Transform data to Eigenvector space
function eigOrient(self,~,~)
if self.img.check
    self.img.eigOrient;
    self.dispUDslice;
    self.dispUDhist;
end