% CMIclass function 
function setImg(self,iimg,labels,fov)
% Set image matrix manually
% Input: iimg = 3D or 4D image matrix

% First set ImageClass matrix:
[d(1),d(2),d(3),d(4)] = size(iimg);
self.img.setMat(iimg,labels,fov./d(1:3));

% Not update image-related CMIclass properties:
self.prmcheck = false;
self.vec = 1;
self.bgvec = 1;
self.orient = 3;
self.slc = round(self.img.dims(1:3)/2);
[vmin,vmax] = self.img.getColorMinMax;
self.clim = [vmin vmax];

% Now update the CMI GUI:
self.GUIupdate;
