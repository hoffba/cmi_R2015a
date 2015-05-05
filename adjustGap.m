% Script to adjust gapped data, inserting blank slices
function adjustGap(cmiObj,targetdz)

ovoxsz = cmiObj.img.voxsz;
dnew = floor(ovoxsz(3)/targetdz);
% vnew = ovoxsz(3)/dnew;
d = cmiObj.img.dims;
fov = d(1:3).*ovoxsz;
d(3) = dnew * d(3);
ind = dnew/2:dnew:d(3);

omask = cmiObj.img.mask.mat;

% Set new image:
iimg = zeros(d)-2000;
iimg(:,:,ind,:) = cmiObj.img.mat;
cmiObj.setImg(iimg,cmiObj.img.labels,fov);

% Set new mask:
if ~isempty(omask)
    imask = false(d(1:3));
    imask(:,:,ind) = omask;
    cmiObj.img.mask.merge('replace',imask);
end
