% ImageClass function
% Reconstruct CT images
function reconCT(self,opts)
if self.check
    answer = questdlg('Is a sinogram loaded?','Sinogram Check','Yes');
    if nargin==1
        strs = {};
        defs = {};
        answer = inputdlg(strs,'CT Reconstruction',1,defs);
    end
    
end