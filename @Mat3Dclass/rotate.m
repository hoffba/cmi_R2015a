% Mat3Dclass function
% Rotate 2D Image (CLOCKWISE)
function rotate(self,dim,a)
% Rotates in current dims([1 2]) without cropping or changing voxel dimensions
% ** Assumes equal voxel spacing within slice plane (otherwise may distort the image)
if (nargin == 3) && self.check && ~isempty(a) && isnumeric(a) && (a~=0) && ismember(dim,1:3)
    
    if isa(self,'ImageClass')
        interpm = 'cubic';
    else
        interpm = 'nearest';
    end
    
    % Determine rotational dimension order:
    switch dim
        case 1 % Row
            p = [1,3,2,4];
        case 2 % Col
            p = [3,1,2,4];
        case 3 % Slc
            p = [2,1,3,4];
    end
    
    % Create transform matrix:
    fov = (self.dims(1:3)-1).*self.voxsz/2;
    R = imref3d(self.dims(1:3),fov(2)*[-1,1],fov(1)*[-1,1],fov(3)*[-1,1]);
    A = [ cosd(a) , -sind(a) , 0 , 0 ;...
          sind(a) , cosd(a)  , 0 , 0 ;...
          0       , 0        , 1 , 0 ;...
          0       , 0        , 0 , 1 ];
    tform = affine3d(A(p,p));
    
    % Loop over 4D:
    nvec = self.dims(4);
    h = waitbar(0,sprintf('Rotating each image: 0 (/%u)',nvec));
    for ivec = 1:nvec
        self.mat(:,:,:,ivec) = imwarp(self.mat(:,:,:,ivec),R,tform,interpm,'OutputView',R);
        waitbar(ivec/nvec,h,sprintf('Rotating each image: %u (/%u)',ivec,nvec));
    end
    close(h)
    
    
    if isprop(self,'mask') && isa(self.mask,'MaskClass') && isvalid(self.mask)
        self.mask.rotate(dim,a);
    end
    if isprop(self,'prm') && isa(self.prm,'PRMclass') && isvalid(self.prm)
        self.prm.rotate(dim,a);
    end
end