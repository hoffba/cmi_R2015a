% CMIclass function
% Allows user to manually align images within a series to first image
function manreg(self)
if self.img.check && self.img.dims(4)>1
    baseimg = self.img.mat(:,:,self.slc(self.orient),1);
    origimg = self.img.mat(:,:,self.slc(self.orient),self.vec);
    dims = size(baseimg);
    
end