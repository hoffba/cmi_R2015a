% Function: coilCombine
%
% Description: combine multi-coil image sequences
%
% Based on: Walsh DO, Gmitro AF, Marcellin MW. Adaptive reconstruction of
% phased array MR imagery. Magn Reson Med 2000;43:682-690
%
% Parameters:
% im1: the multi-coil images (size [nx,ny,nz,ncoils,nTE])
%
% Returns:
% im2: the coil-combined images (size [nx,ny,nz,1,nTE])
%
% Author: Diego Hernando
% Date created: August 13, 2011
% Date last modified: December 8, 2011

function im2 = coilCombine( im1 )

[sx,sy,sz,C,N] = size(im1);

% Let's make the coil dimension last
im1 = permute(im1,[1 2 3 5 4]);

% Set filter size
filtsize = 7;
fkern = ones(([sx,sy,sz]>filtsize)*(filtsize-1) + 1);

% Initialize
im2 = zeros(sx,sy,sz,1,N);
Rs = zeros(sx,sy,sz,C,C);

% Get correlation matrices
for kc1=1:C
    for kc2=1:C
        for kn=1:N
            Rs(:,:,:,kc1,kc2) = Rs(:,:,:,kc1,kc2) + convn(im1(:,:,:,kn,kc1).*conj(im1(:,:,:,kn,kc2)),fkern,'same');
        end
    end
end

% Compute and apply filter at each voxel
for kx = 1:sx
    for ky = 1:sy
        for kz = 1:sz
            % $$$     [U,S] = eig(squeeze(Rs(kx,ky,:,:)));
            % $$$     s = diag(S);
            % $$$     [maxval,maxind] = max(abs(s));
            % $$$     myfilt = U(:,maxind);
            % $$$     im2(kx,ky,:) = myfilt'*reshape(squeeze(im1(kx,ky,:,:)).',[C N]);

            % Change suggested by Mark Bydder
            [U,~] = svd(squeeze(Rs(kx,ky,kz,:,:)));
            myfilt = U(:,1);
            im2(kx,ky,kz,1,:) = myfilt'*reshape(squeeze(im1(kx,ky,kz,:,:)).',[C N]);
        end
    end
end

% In case the input data are single
if isa(im1,'single')
    im2 = single(im2);
end
