% CMIclass function
% Performs perfusion MRI calculation
function calcFit(self,~,~)
    if self.img.check && self.img.mask.check
        self.img.calcFit;
        % Update image values
        self.clim = [squeeze(min(min(min(self.img.mat)))),...
                     squeeze(max(max(max(self.img.mat))))];
        % Now update GUI objects with new image
        nv = self.img.dims(4);
        set(self.h.slider_4D,'Value',1,'SliderStep',[1,1]/(nv-1),'Max',nv);
        set(self.h.text_nv,'String',num2str(nv));
        self.setVec(1);
    end
end

