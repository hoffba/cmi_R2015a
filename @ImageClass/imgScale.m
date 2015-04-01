% ImageClass function
% Scale current image values (y = m*x + b)
function imgScale(self,vec,m,b)
if self.check
    if (nargin == 4)
        if isnumeric(m) && isnumeric(b) && ~any(vec<1) && ~any(vec>self.dims(4))
            %m = m(1); b = b(1); % only use one value
            vec = uint8(vec); % make sure they're ints
            for i = 1:length(vec)
                tm = self.scaleM(vec(i));
                tb = self.scaleB(vec(i));
                self.mat(:,:,:,vec(i)) = m(i)/tm * (self.mat(:,:,:,vec(i)) - tb) + b(i);
                self.valExt(vec(i),:) = m(i)/tm * (self.valExt(vec(i),:) - tb) + b(i);
                self.scaleM(vec(i)) = m(i);
                self.scaleB(vec(i)) = b(i);
                self.thresh(vec(i),:) = (self.thresh(vec(i),:) - tb) * m(i)/tm + b(i);
            end
        end
    end
end