% ImageClass function
function volCorrectMRI(self,vec)

self.mat(:,:,:,vec) = vol_homocor_bah_3D(self.mat(:,:,:,vec));