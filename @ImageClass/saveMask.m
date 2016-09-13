% ImageClass function
% Save mask
function status = saveMask(self,fname)
if (nargin<2) || isempty(fname)
    fname = fullfile(self.dir,[self.name,'_VOI.mhd']);
end
status = self.mask.save(self.voxsz,fname);