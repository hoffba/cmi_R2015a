% ImageClass function
function [X,Y,Z] = getImageCoords(self)

n = cross(self.dircos(1:3),self.dircos(4:6));
M = [reshape(self.dircos,3,2),n',self.slcpos';zeros(1,3),1];

d = self.dims(1:3);
[X,Y,Z] = meshgrid((0:d(2)-1)*self.voxsz(2),(0:d(1)-1)*self.voxsz(1),(0:d(3)-1)*self.voxsz(3));
I = M * [reshape(cat(4,X,Y,Z),prod(d),3)';ones(1,prod(d))];
X = reshape(I(1,:),d);
Y = reshape(I(2,:),d);
Z = reshape(I(3,:),d);

% d = self.dims(1:3);
% v = self.voxsz;
% vd = (d-1)/2 .* v;
% [X,Y,Z] = meshgrid(-vd(2):v(2):vd(2),-vd(1):v(1):vd(1),-vd(3):v(3):vd(3));
