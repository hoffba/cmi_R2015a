% ImageClass function
function vol = getVol(self)
vol = sum(self.mask.mat)*prod(self.voxsz);