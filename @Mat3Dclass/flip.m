% Mat3Dclass function
% Flip Image in true dimension "dim"
function flip(self,dim)
if (nargin == 2) && ~ischar(dim) && ~isempty(dim) && ~any((dim > 3) | (dim < 1))
    for i = 1:length(dim)
        if self.check
            self.mat = flip(self.mat,dim(i));
        end
    end
    if isprop(self,'mask') && isa(self.mask,'MaskClass')
        self.mask.flip(dim);
    end
    if isprop(self,'prm') && isa(self.prm,'PRMclass')
        self.prm.flip(dim);
    end
end