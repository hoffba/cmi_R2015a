% CMIclass function
function setDCEopts(self,~,~)
if isa(self.img.model,'DCEclass')
    self.img.model.setDCEopts;
end