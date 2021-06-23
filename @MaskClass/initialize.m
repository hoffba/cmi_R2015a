% MaskClass function
% Re-initialize the mask
function initialize(self,dims)
if (nargin==2) && isnumeric(dims) && (length(dims)==3)
    self.mat = false(dims);
    self.dims = [dims,1];
    self.check = false;
else
    self.mat = [];
    self.dims = zeros(1,4);
    self.check = false;
end