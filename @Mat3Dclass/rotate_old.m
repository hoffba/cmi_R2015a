% Mat3Dclass function
% Rotate 2D Image (CLOCKWISE)
function rotate(self,dim,a)
% Rotates in current dims([1 2]) without cropping or changing voxel dimensions
% ** Assumes equal voxel spacing within slice plane (otherwise may distort the image)
if (nargin == 3) && self.check && ~isempty(a) && isnumeric(a) && (a~=0) ...
        && ismember(dim,1:3)
    
    % Determine rotational dimension order:
    switch dim
        case 1 % Row
            p = [2,3,1];
        case 2 % Col
            p = [1,3,2];
        case 3 % Slc
            p = [1,2,3];
    end
    
    % Create transform matrix:
    A = [ cosd(a) , -sind(a) , 0 ;...
          sind(a) , cosd(a)  , 0 ;...
          0       , 0        , 1 ];
    T = maketform('affine',A);
    
    % Gather image to rotate
    tmat = permute(self.mat,[p,4]);
    [tdims(1),tdims(2),tdims(3),tdims(4)] = size(tmat);
    pxsz = self.voxsz(p(1:2));
    ext = (tdims(1:2) - 1)/2 .* pxsz;
    
    % Start rotating each slice:
    ntot = prod(tdims(3:4));
    h = waitbar(0,['Rotating each slice: 0 (/' num2str(ntot) ')']);
    count = 0;
    for islc = 1:tdims(3)
        for ivec = 1:tdims(4)
            timg = imtransform(tmat(:,:,islc,ivec),T,'UData',ext(2)*[-1,1],'VData',ext(1)*[-1,1]);
            if count==0 % Performed on first pass to initialize the new image matrix
                d = self.dims;
                d(p(1:2)) = size(timg);
                self.mat = zeros(d);
                self.dims = d;
%                 self.setMat(zeros(d));
            end
            self.setSlice(timg,dim,ivec,islc);
            count = count + 1;
            waitbar(count/ntot,h,['Rotating each slice: ' num2str(count)...
                ' (/' num2str(ntot) ')']);
        end
    end
    close(h)
    
    
    if isprop(self,'mask') && isa(self.mask,'MaskClass') && isvalid(self.mask)
        self.mask.rotate(dim,a);
    end
    if isprop(self,'prm') && isa(self.prm,'PRMclass') && isvalid(self.prm)
        self.prm.rotate(dim,a);
    end
end