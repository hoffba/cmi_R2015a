% MaskClass function
function setMat(self,imat)
% Manually set the matrix
[d(1),d(2),d(3)] = size(imat);
if all(d==self.dims(1:3))
    self.mat = imat;
    self.check = true;
end