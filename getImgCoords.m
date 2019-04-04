function [x,y,z] = getImgCoords(d,voxsz,dircos,offset)
% Returns meshgrid of image voxel coordinates

np = prod(d);
M = [ reshape(dircos,3,3),offset(:) ; zeros(1,3),1];
[y,x,z] = meshgrid((0:d(2)-1)*voxsz(2),...
                   (0:d(1)-1)*voxsz(1),...
                   (0:d(3)-1)*voxsz(3));
xyz = (M * [x(:),y(:),z(:),ones(np,1)]')';
x = reshape(xyz(:,1),d);
y = reshape(xyz(:,2),d);
z = reshape(xyz(:,3),d);