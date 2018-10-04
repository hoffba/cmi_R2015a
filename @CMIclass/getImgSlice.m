% CMIclass function
% Retrieve current slice
function islc = getImgSlice(self,chk,vec)
% Inputs: chk = 'mask', 'img', or 'prm'
%         vec = (optional) image 4D index
switch chk
    case 'mask'
        islc = self.img.mask.getSlice(self.orient,1,self.slc(self.orient));
    case 'img'
        if nargin<3
            vec = self.vec;
        elseif (vec<1) || (vec>self.img.dims(4))
            error(['Requested vec#: ',num2str(vec),' - Out of bounds.']);
        end
        islc = self.img.getSlice(self.orient,vec,self.slc(self.orient));
    case 'prm'
        islc = self.img.prm.getSlice(self.orient,1,self.slc(self.orient));
    otherwise
        islc = [];
end
