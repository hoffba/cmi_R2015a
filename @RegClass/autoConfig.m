% RegClass function
function autoConfig(self, ~, ~)
% Function for automatically determining SP_a parameter for warping
%   based on reference image dimensions

if self.cmiObj(1).img.check
    
    Tstr = self.elxObj.getPar(self.ind,'Transform');
    if strcmp(Tstr,'BSplineTransform')
        [A,alph,gam,g,n] = self.elxObj.getPar(self.ind,'SP_A','SP_alpha',...
            'FinalGridSpacingInVoxels','GridSpacingSchedule','NumberOfResolutions');
        v = self.cmiObj(1).img.voxsz*1000;
        a = round(min(repmat((v.*gam)',1,n).*reshape(g,[3,n])) .* A.^alph / 0.4);
    else
        a = self.elxObj.getPar(self.ind,'SP_a');
    end
    self.setElxPar(self.ind,'SP_a',a);
end