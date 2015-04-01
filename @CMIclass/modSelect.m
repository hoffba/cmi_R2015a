% CMIclass function
% Select which model to perform fit with
function modSelect(self,~,~)

self.img.model.setModel;
set(self.h.analysis_modsel,'Label',...
    ['Model: ',self.img.model.getModName]);