% ImageClass function
% Set image voxel dimensions
function setVoxSz(self,sz)
if (nargin == 2) && self.check && (numel(sz)==3)
    self.voxsz = sz(:)';
end