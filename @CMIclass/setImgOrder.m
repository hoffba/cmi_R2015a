% CMIclass function
function setImgOrder(self,x,~)

if isa(x,'')
elseif isnumeric(x)
    x = x(:);
end
self.img.setImgOrder(x);
self.dispUDsice;