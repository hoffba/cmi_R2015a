% CMIclass function
function calcLinLSF(self,~,~)
%calcLinLSF Summary of this function goes here
%   Detailed explanation goes here

self.img.expFit;
self.bgvec = 1;
[vmin,vmax] = self.img.getColorMinMax;
self.clim = [vmin,vmax];
nv = self.img.dims(4);
set(self.h.text_nv,'String',num2str(nv));
set(self.h.slider_4D,'Value',1,'Max',nv,'SliderStep',[1,1]/(nv-1));
self.setVec(1);

