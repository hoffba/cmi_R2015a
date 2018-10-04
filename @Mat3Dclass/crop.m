% Mat3Dclass function
% Crop Image to current zoom
function crop(self,dim,lims)
% Inputs:
%   dim = [n x 1] column vector of dimensions to crop (1, 2, or 3)
%   lims = [n x 2] matrix with rows [min max] describing matrix crop limits

if (nargin == 3) && all(ismember(dim,1:3)) ...
        && isnumeric(lims) && (size(dim,1)==size(lims,1)) && (size(lims,2)==2)
    lims = round(lims); % Make sure they're integer values
    if ~isempty(self.mat)
        % Crop each dimension
        iy = 1:self.dims(1);
        ix = 1:self.dims(2);
        iz = 1:self.dims(3);
        for i = 1:length(dim)
            switch dim(i)
                case 1
                    iy((iy<lims(i,1))|(iy>lims(i,2))) = [];
                case 2
                    ix((ix<lims(i,1))|(ix>lims(i,2))) = [];
                case 3
                    iz((iz<lims(i,1))|(iz>lims(i,2))) = [];
            end
        end
        self.mat = self.mat(iy,ix,iz,:);
        self.slcpos = self.getPosition('Coordinates',[iy(1),ix(1),iz(1)]);
    end
    self.dims(dim) = diff(lims,1,2)+1;
    if isprop(self,'mask') && isa(self.mask,'MaskClass') && isvalid(self.mask)
        self.mask.crop(dim,lims);
    end
    if isprop(self,'prm') && isa(self.prm,'PRMclass') && isvalid(self.prm)
        self.prm.crop(dim,lims);
    end
end