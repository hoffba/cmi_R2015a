function [img,label,fov,info] = readNIFTI(varargin)
%% reads Nifti *.nii file
fname = varargin{1};

[~, Nfile, ~] = fileparts(fname);

img = double(niftiread(fname));
info = niftiinfo(fname);

img = rot90(fliplr(img));

fov = info.ImageSize.*info.PixelDimensions;
label = {Nfile};

