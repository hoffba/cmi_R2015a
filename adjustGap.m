% Script to adjust gapped data, inserting blank slices
function adjustGap(cmiObj,targetdz,fillval)

if nargin<3
    fillval = 0;
end

dnew = floor(cmiObj.img.voxsz(3)/targetdz);
d = cmiObj.img.dims;
fov = d(1:3).*cmiObj.img.voxsz;
d(3) = dnew * d(3);
ind = dnew/2:dnew:d(3);

omask = cmiObj.img.mask.mat;

% Set new image:
iimg = zeros(d)+fillval;
iimg(:,:,ind,:) = cmiObj.img.mat;
cmiObj.setImg(iimg,cmiObj.img.labels,fov);

% Set new mask:
if ~isempty(omask)
    imask = false(d(1:3));
    imask(:,:,ind) = omask;
    cmiObj.img.mask.merge('replace',imask);
end
