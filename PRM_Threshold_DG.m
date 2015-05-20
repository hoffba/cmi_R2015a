data=cmiObj0.img.mat;
mask=cmiObj0.img.mask.mat;
wk0=data(:,:,:,1);wk10=data(:,:,:,2);
vwk0=wk0(mask);vwk10=wk10(mask);
PRM_threshold2([vwk0 vwk10])