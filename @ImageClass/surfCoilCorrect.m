% ImageClass function
function surfCoilCorrect(self,vec,thresh)
% Corrects signal dropoff in MRI images

if nargin<3 || ~isnumeric(thresh)
    thresh = str2double(inputdlg('Input Threshold:','Threshold',1,{'0.02'}));
end
self.mat(:,:,:,vec) = vol_homocor_bah_3D(self.mat(:,:,:,vec),thresh);
