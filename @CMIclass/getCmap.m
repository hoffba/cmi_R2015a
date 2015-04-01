% CMIcalss function
% Retrieve current colormap
function map = getCmap(self)
if self.bgcheck
    map = self.bgcmap;
else
    map = self.cmap;
end