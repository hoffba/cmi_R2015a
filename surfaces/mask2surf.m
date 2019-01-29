function [fv,x,y,z] = mask2surf(BW,voxsz,V)

d = size(BW);
if nargin<2
    voxsz = ones(1,3);
end
if nargin<5
    V = {};
else
    V = {V};
end

% np = prod(d);
[y,x,z] = meshgrid((0:d(2)-1)*voxsz(2),...
                   (0:d(1)-1)*voxsz(1),...
                   (0:d(3)-1)*voxsz(3));
% xyz = M * [reshape(cat(4,x,y,z),np,3)';ones(1,np)];
% xyz = reshape(xyz(1:3,:)',[d,3]);
% x = xyz(:,:,:,1);
% y = xyz(:,:,:,2);
% z = xyz(:,:,:,3);
% clear xyz

fv = isosurface(x,y,z,smooth3(double(BW)*1000,'gaussian',5,1.5),500,V{:});
fv = smoothpatch(fv,1,5,1,[]);
fv = reducepatch(fv);
