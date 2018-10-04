% MaskClass function
% Re-initialize the mask
function initialize(self,dims,voxsz,dircos,slcpos)
if nargin<2 || (numel(dims)~=3) || ~isnumeric(dims)
    dims = zeros(1,3);
end
if nargin<3 || (numel(voxsz)~=3) || ~isnumeric(voxsz)
    voxsz = ones(1,3);
end
if nargin<4 || (numel(dircos)~=6) || ~isnumeric(dircos)
    dircos = [1 0 0 0 1 0];
end
if nargin<5 || (numel(slcpos)~=3) || ~isnumeric(slcpos)
    slcpos = zeros(1,3);
end
self.mat = [];
self.dims = [dims(:)',1];
self.voxsz = voxsz;
self.dircos = dircos;
self.slcpos = slcpos;
self.check = false;

