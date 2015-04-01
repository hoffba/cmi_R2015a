% CMIclass function
function CTdensityCal(self,opt,~)
if (nargin==1) || ishandle(opt)
    opt = [];
end
[m,b] = self.img.CTphantCal(self.vec,opt);
self.clim(self.vec,:) = self.clim(self.vec,:) * m + b;
set(self.h.edit_cmin,'String',num2str(self.clim(self.vec,1)))
set(self.h.edit_cmax,'String',num2str(self.clim(self.vec,2)))
self.dispUDslice;