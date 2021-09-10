% RegClass function
function swapCMIdata(self,keepname)
% Swaps data in CMIclass objects
%  * useful for registering serial images
if nargin==1
    keepname = false;
end

dimchk = self.cmiObj(1).img.dims(4)==self.cmiObj(2).img.dims(4);

% Grab temporary fixed image values:
timg = self.cmiObj(1).img.mat;
tfov = self.cmiObj(1).img.voxsz.*self.cmiObj(1).img.dims(1:3);
tmask = self.cmiObj(1).img.mask.mat;
torient = self.cmiObj(1).img.orient;

% Optional name info swap:
tlab = { self.cmiObj(1).img.labels , self.cmiObj(2).img.labels };
tname = { self.cmiObj(1).img.name , self.cmiObj(2).img.name };
tdir = { self.cmiObj(1).img.dir , self.cmiObj(2).img.dir };

if ~keepname
    tname = tname([2,1]);
    tdir = tdir([2,1]);
    if dimchk
        tlab = tlab([2,1]);
    end
end

% Set new fixed image to old moving image values:
self.cmiObj(1).setImg(self.cmiObj(2).img.mat,...
                      tlab{1},...
                      self.cmiObj(2).img.voxsz.*self.cmiObj(2).img.dims(1:3),...
                      self.cmiObj(2).img.orient,...
                      fullfile(tdir{1},tname{1}));
self.cmiObj(1).img.mask.merge('replace',self.cmiObj(2).img.mask.mat);
self.cmiObj(1).setSlice(self.cmiObj(1).slc(self.cmiObj(1).orient));

% Set new moving image to old fixed image values:
self.cmiObj(2).setImg(timg,tlab{2},tfov,fullfile(tdir{2},torient,tname{2}));
self.cmiObj(2).img.mask.merge('replace',tmask);
self.cmiObj(2).setSlice(self.cmiObj(2).slc(self.cmiObj(1).orient));