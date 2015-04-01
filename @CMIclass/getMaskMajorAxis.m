% CMIclass function
function D = getMaskMajorAxis(self,~,~)

D = self.img.getMaskMajorAxis(self.orient);
if nargout==0
    disp(['Major Axis Diameter = ',num2str(D)])
end