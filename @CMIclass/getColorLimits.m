% CMIcalss function
% Retrieve color limits for current image
function [vmin,vmax] = getColorLimits(self)
if self.bgcheck
    tvec = self.bgvec;
else
    tvec = self.vec;
end
vmin = self.clim(tvec,1);
vmax = self.clim(tvec,2);