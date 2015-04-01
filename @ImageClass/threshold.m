% ImageClass function
% Threshold Image
function threshold(self,vec,val)
% Inputs: vec = [1xn] vector of image indices
%         val = [2xn] matrix of [min;max] values for each vector
if (nargin == 3) && isnumeric(val) && isnumeric(vec)
    vec = round(vec(:));
    nv = length(vec);
    if (size(val,1) == 2) && all(vec>0) && all(vec<=self.dims(4))...
            && (nv == size(val,2)) && all(val(1,:)<val(2,:))
        for i = 1:nv
            self.thresh(vec(i),:) = val(:,i);
        end
    end
end