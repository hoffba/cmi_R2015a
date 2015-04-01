% ImageClass function
% Save mask
function status = saveMask(self,fname)
if (nargin<2) || isempty(fname)
    [~,fname,~] = fileparts(self.name);
    fname = [fname '_VOI'];
end
status = self.mask.save(self.voxsz,fname);