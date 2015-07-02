% CMIcalss function
% Retrieve current slice
function islc = getImgSlice(self,chk)
% Inputs: chk = 'mask', 'img', or 'prm'
switch chk
    case 'mask'
        islc = self.img.mask.getSlice(self.orient,1,self.slc(self.orient));
    case 'img'
        islc = self.img.getSlice(self.orient,self.vec,self.slc(self.orient));
    case 'prm'
        islc = self.img.prm.getSlice(self.orient,1,self.slc(self.orient));
    otherwise
        islc = [];
end
