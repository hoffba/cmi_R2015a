% CMIclass function
function update4Dslider(self)

nv = self.img.dims(4);
set(self.h.slider_4D,'Max',nv,'SliderStep',[1,1]/(nv-1))