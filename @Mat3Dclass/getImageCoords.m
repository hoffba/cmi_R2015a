% Mat3Dclass function
function [X,Y,Z] = getImageCoords(self)
d = self.dims(1:3);
v = self.voxsz;
vd = (d-1)/2 .* v;
[X,Y,Z] = meshgrid(-vd(2):v(2):vd(2),-vd(1):v(1):vd(1),-vd(3):v(3):vd(3));
