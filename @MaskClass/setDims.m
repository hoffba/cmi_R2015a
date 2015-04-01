% MaskClass function
% Set desired mask dimensions
function setDims(self,d)
if (nargin==2) && isnumeric(d) && (length(d)==3)
    self.dims = [d 1];
end