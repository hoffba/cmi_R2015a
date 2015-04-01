% Mat3Dclass
function permuteMat(self,order)
% Allows for permutation of first 3 dimensions
% Input = some permutation of [1,2,3]
%   - no other numbers!
%   - includes all three!
if (nargin == 2) && (length(order)==3) && self.check
    order = [round(order),4];
    if all(sort(order)==(1:4))
        self.mat = permute(self.mat,order);
        self.dims = self.dims(order);
    end
    if isa(self,'ImageClass')
        self.setVoxSz(self.voxsz(order(1:3)));
    end
    if isprop(self,'mask') && isa(self.mask,'MaskClass') && isvalid(self.mask)
        self.mask.permuteMat(order(1:3));
    end
    if isprop(self,'prm') && isa(self.prm,'PRMclass') && isvalid(self.prm)
        self.prm.permuteMat(order(1:3));
    end
end