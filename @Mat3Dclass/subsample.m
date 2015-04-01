% Mat3Dclass function
% Resample image/VOI matrix
function stat = subsample(self,ix,iy,iz)
% Inputs:   ix = dim1 indices to keep (empty keeps all)
%           iy = dim2 
%           iz = dim3 
%   Warning!!! Subsampling keeps image voxel size, so true spatial 
%                   coordinates are not maintained

stat = false;
if self.check && (nargin==4)
    if isempty(ix)
        ix = 1:self.dims(1);
    end
    if isempty(iy)
        iy = 1:self.dims(1);
    end
    if isempty(iz)
        iz = 1:self.dims(1);
    end
    self.mat = self.mat(ix,iy,iz,:);
    self.dims(1) = length(ix);
    self.dims(2) = length(iy);
    self.dims(3) = length(iz);
    stat = true;
    if isa(self,'ImageClass') && self.mask.check
        if self.mask.check
            self.mask.subsample(ix,iy,iz);
        end
        if self.prm.check
            self.prm.subsample(ix,iy,iz);
        end
    end
end

