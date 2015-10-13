% ImageClass function
function shiftMat(self,v,dx)
% Manually shift the matrix

if self.check && all(ismember(v,1:self.dims(4))) && (length(dx)==3)
    self.mat(:,:,:,v) = circshift(self.mat(:,:,:,v),round(dx));
end