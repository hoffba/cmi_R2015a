% Mat3Dclass function
% Apply affine spatial transform to image
function affTxF(self,T,interpm)
% interpm = 'linear', 'nearest'
% T       = affine transformation matrix
if self.check && (nargin>1)
    % ~~~~~ Start transformation ~~~~~~~~~~~~~~~~~
    tic
    d = self.dims;
    n = prod(d(1:3));
    for i = 1:d(4)
        timg = self.mat(:,:,:,i);
        if isa(self,'MaskClass')
            timg = double(timg);%getOrigData(i);
        end
        tmin = min(timg(:));
        if (nargin<3) || isempty(interpm) || ~any(strcmp(interpm,{'linear','nearest'}))
            interpm = 'linear';
        end
        [X,Y,Z] = meshgrid(0:d(2)-1,0:d(1)-1,0:d(3)-1);
        C=[X(:)';Y(:)';Z(:)';ones(1,n)];
        Ci=T\C;
        Xi=reshape(Ci(1,:),d(1:3));
        Yi=reshape(Ci(2,:),d(1:3));
        Zi=reshape(Ci(3,:),d(1:3));
        timg = interp3(timg,Xi,Yi,Zi,interpm);
        timg(isnan(timg)) = tmin;
        if isa(self,'MaskClass')
            timg = logical(timg);
        end
        self.mat(:,:,:,i) = timg;
    end
    toc
    % ~~~~~ End transformation ~~~~~~~~~~~~~~~~~
end