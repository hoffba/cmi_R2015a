% CMIcalss function
% Retrieve intensity min/max for current image
function [vmin,vmax] = getColorMinMax(self)
if self.bgcheck
    tvec = self.bgvec;
else
    tvec = self.vec;
end
[vmin,vmax] = self.img.getColorMinMax(tvec);