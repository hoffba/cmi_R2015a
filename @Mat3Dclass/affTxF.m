% Mat3Dclass function
% Apply affine spatial transform to image
function affTxF(self,T,interpm)
% interpm = 'linear', 'nearest'
% T       = affine transformation matrix
if self.check && (nargin>1)
    maskchk = isa(self,'MaskClass');
    if (nargin<3) || isempty(interpm)
        interpm = 'linear';
    end
    if maskchk
        interpm = 'nearest';
    end
    % Image coordinates:
    
    % ~~~~~ Start transformation ~~~~~~~~~~~~~~~~~
    tic
    d = self.dims;
    n = prod(d(1:3));
    for i = 1:d(4)
        timg = self.mat(:,:,:,i);
        if maskchk
            timg = double(timg);
        end
        tmin = min(timg(:));
        [X,Y,Z] = self.getImageCoords;
        C=[X(:)';Y(:)';Z(:)';ones(1,n)];
%         Ci=T\C;
        Ci = (C' * T)';
        Xi=reshape(Ci(1,:),d(1:3));
        Yi=reshape(Ci(2,:),d(1:3));
        Zi=reshape(Ci(3,:),d(1:3));
        timg = interp3(X,Y,Z,timg,Xi,Yi,Zi,interpm);
        timg(isnan(timg)) = tmin;
        if maskchk
            timg = logical(timg);
        end
        self.mat(:,:,:,i) = timg;
    end
    toc
    % ~~~~~ End transformation ~~~~~~~~~~~~~~~~~
end