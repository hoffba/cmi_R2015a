% ImageClass function
% Load mask
function status = loadMask(self,fname,optsel)
if nargin < 2
    fname = [];
end
if nargin < 3
    optsel = '';
end

if self.check % Only load a mask if image exists
    status = self.mask.load(self.dims(1:3),self.voxsz,fname,optsel);
else
    status = false;
end