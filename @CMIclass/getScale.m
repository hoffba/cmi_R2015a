% CMIcalss function
% Retrieve current image scaling values (m,b)
function [m,b] = getScale(self)
m = self.img.scaleM(self.vec);
b = self.img.scaleB(self.vec);