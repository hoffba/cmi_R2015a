% Mat3Dclass function
function img = getSlice(self,vdim,vec,slc)
% Retrieve slice from 3D image matrix

img = [];
if (nargin==4) && ismember(vdim,1:3) ...
        && ismember(vec,1:self.dims(4)) ...
        && ismember(slc,1:self.dims(vdim))
    d = self.dims(1:3); d(vdim) = [];
    img = zeros(d);
    if self.check
        switch vdim
            case 1 % ROW
                img(:) = self.mat(slc,:,:,vec);
            case 2 % COL
                img(:) = self.mat(:,slc,:,vec);
            case 3 % SLC
                img(:) = self.mat(:,:,slc,vec);
        end
    end
end