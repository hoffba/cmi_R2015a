% CMIclass function
% Retrieve slice limits
function nslices = getNSlices(self)
nslices = self.img.dims(self.orient);