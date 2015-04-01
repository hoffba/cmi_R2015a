% Mat3Dclass function
function setSlice(self,tmat,vdim,vec,slc)
% Set slice into 3D image matrix

if (nargin==5) && ismember(vdim,1:3) ...
        && ismember(vec,1:self.dims(4)) ...
        && ismember(slc,1:self.dims(vdim))
    d = self.dims(1:3); d(vdim) = [];
    if all(size(tmat)==d)
        if ~self.check
            self.mat = false(self.dims);
        end
        switch vdim
            case 1 % ROW
                self.mat(slc,:,:,vec) = tmat;
            case 2 % COL
                self.mat(:,slc,:,vec) = tmat;
            case 3 % SLC
                self.mat(:,:,slc,vec) = tmat;
        end
        if any(self.mat(:))
            self.check = true;
        else
            self.check = false;
            self.mat = [];
        end
    end
end