% ImageClass function
% Load mask
function status = loadMask(self,fname)
if (nargin == 1)
    fname = [];
end
if self.check % Only load a mask if image exists
    status = self.mask.load(self.dims(1:3),self.voxsz,fname);
else
    status = false;
end