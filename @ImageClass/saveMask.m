% ImageClass function
% Save mask
function status = saveMask(self,fname)
if (nargin<2) || isempty(fname)
    fname = [self.name,'_VOI.nii.gz'];
end
status = self.mask.save(self.voxsz,self.orient,fname);