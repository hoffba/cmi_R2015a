% CMIclass function
function imgMath(self,~,~,~)
% Performs calculation on current image
if self.img.check
    ok = self.img.imgMath;
    if ok
        nvec = size(self.clim,1)+1 : self.img.dims(4);
        self.clim(nvec,:) = [squeeze(min(min(min(self.img.mat(:,:,:,nvec))))),...
                             squeeze(max(max(max(self.img.mat(:,:,:,nvec)))))];
        d = self.img.dims(4);
        set(self.h.text_nv,'String',num2str(d));
        set(self.h.slider_4D,'Max',d,'SliderStep',[1,1]/(d-1));
        self.setVec(d);
    end
end