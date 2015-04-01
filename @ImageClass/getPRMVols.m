% ImageClass function
function vols = getPRMVols(self)
if self.prm.check
    vols = zeros(self.prm.nprm,1);
    for i = 1:self.prm.nprm
        vols(i) = sum(self.prm.mat(:)==i);
    end
    vols = vols * prod(self.voxsz);
end
