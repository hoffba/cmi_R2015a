% CMIclass function
% Transform data to Eigenvector space
function eigOrient(self,x,~)
if self.img.check
    if nargin==1 || isa(x,'matlab.ui.container.Menu')
        x = [];
    end
    self.img.eigOrient(x);
    self.dispUDslice;
    self.dispUDhist;
end