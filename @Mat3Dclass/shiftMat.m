% Mat3Dclass function
function shiftMat(self,v,dx)
% Shift (circular) the matrix

if self.check && all(ismember(v,1:self.dims(4))) && (length(dx)==3)
    
    % Shift the image matrix
    self.mat(:,:,:,v) = circshift(self.mat(:,:,:,v),round(dx));
    
    if isprop(self,'mask')
        self.mask.shiftMat(1,dx);
    end
end