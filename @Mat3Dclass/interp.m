% Mat3Dclass function
% Interpolate Image
function interp(self,nd)
if (nargin == 2) && (length(nd)==3)
    if self.check
        % Make sure all 3 dimensions are included
        nd = nd(:)';
        d = self.dims(1:3);
        dd = (d(:)'./nd);
        [xi,yi,zi] = meshgrid((dd(2) : dd(2) : d(2))+0.5-dd(2)/2,...
                              (dd(1) : dd(1) : d(1))+0.5-dd(1)/2,...
                              (dd(3) : dd(3) : d(3))+0.5-dd(3)/2);
        lchk = islogical(self.mat);
        if lchk
            self.mat = single(self.mat);
        end
        tmat = zeros([nd,self.dims(4)]);
        for i = 1:self.dims(4)
            tmat(:,:,:,i) = interp3(self.mat(:,:,:,i),xi,yi,zi,'linear');
        end
        % Now update the image
        tmat(isnan(tmat)) = 0;
        if lchk
            tmat = logical(round(tmat));
        end
        if isprop(self,'voxsz') && ismethod(self,'setVoxSz')
            self.setVoxSz(self.dims(1:3).*self.voxsz./nd);
        end
        self.dims(1:3) = nd;
        self.mat = tmat;
    end
    self.dims(1:3) = nd;
    if isprop(self,'mask') && isa(self.mask,'MaskClass') && isvalid(self.mask)
        self.mask.interp(nd);
    end
    if isprop(self,'prm') && isa(self.prm,'PRMclass') && isvalid(self.prm)
        self.prm.interp(nd);
    end
end