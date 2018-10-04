function [fv,x,y,z] = mask2surf(BW,voxsz,V)

if nargin<2
    voxsz = ones(1,3);
end
if nargin<3
    V = {};
else
    V = {V};
end

% Center image at origin:
d = size(BW);
lims = (d-1).*voxsz/2;
[y,x,z] = meshgrid(-lims(2):voxsz(2):lims(2),...
                   -lims(1):voxsz(1):lims(1),...
                   -lims(3):voxsz(3):lims(3));

fv = isosurface(x,y,z,smooth3(double(BW)*1000,'gaussian',5,1.5),500,V{:});
fv = smoothpatch(fv,1,5,1,[]);
fv = reducepatch(fv);
