% RegClass function
function swapCMIdata(self)
% Swaps data in CMIclass objects
%  * useful for registering serial images

% Grab temporary fixed image values:
timg = self.cmiObj(1).img.mat;
tfov = self.cmiObj(1).img.voxsz.*self.cmiObj(1).img.dims(1:3);
tlab = self.cmiObj(1).img.labels;
tmask = self.cmiObj(1).img.mask.mat;
tname = self.cmiObj(1).img.name;
tdir = self.cmiObj(1).img.dir;

% Set new fixed image to old moving image values:
self.cmiObj(1).setImg(self.cmiObj(2).img.mat,...
                      self.cmiObj(2).img.labels,...
                      self.cmiObj(2).img.voxsz.*self.cmiObj(2).img.dims(1:3),...
                      fullfile(self.cmiObj(2).img.dir,self.cmiObj(2).img.name));
self.cmiObj(1).img.mask.merge('replace',self.cmiObj(2).img.mask.mat);
self.cmiObj(1).setSlice(self.cmiObj(1).slc(self.cmiObj(1).orient));

% Set new moving image to old fixed image values:
self.cmiObj(2).setImg(timg,tlab,tfov,fullfile(tdir,tname));
self.cmiObj(2).img.mask.merge('replace',tmask);
self.cmiObj(2).setSlice(self.cmiObj(2).slc(self.cmiObj(1).orient));