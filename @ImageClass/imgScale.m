% ImageClass function
% Scale current image values (y = m*x + b)
function imgScale(self,vec,m,b,relchk)
% Inputs:
%   vec = [1xn] vector of 4D image indices to scale
%   m   = [1xn] vector of scale factors
%   b   = [1xn] vector of offsets
%   relchk = (optional) determines whether to scale on the original data (0)
%                       or the currently scaled data (1)
if self.check
    if (nargin >= 4)
        if nargin<5
            relchk = false;
        else relchk = logical(relchk);
        end
        if isnumeric(m) && isnumeric(b) && ~any(vec<1) && ~any(vec>self.dims(4))
            vec = uint8(vec); % make sure they're ints
            for i = 1:length(vec)
                if relchk
                    self.mat(:,:,:,vec(i)) = m(i)*self.mat(:,:,:,vec(i))+b(i);
                    self.valExt(vec(i),:) = m(i)*self.valExt(vec(i),:)+b(i);
                    self.scaleM(vec(i)) = m(i)*self.scaleM(vec(i));
                    self.scaleB(vec(i)) = m(i)*self.scaleB(vec(i))+b(i);
                    self.thresh(vec(i),:) = m(i)*self.thresh(vec(i),:)+b(i);
                else
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
end