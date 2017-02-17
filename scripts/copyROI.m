function copyROI(cmiObj,oslc,slc)


smask = cmiObj.img.mask.getSlice(3,1,oslc);
mask = false(cmiObj.img.dims(1:3));
mask(:,:,slc) = repmat(smask,[1,1,length(slc)]);
cmiObj.img.mask.merge('replace',mask);