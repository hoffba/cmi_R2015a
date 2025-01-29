function status = saveNIFTI(fname,img,label,fov,orient)
% Save image as NIfTI (.nii.gz)

%% Generate NIfTI metadata
d = size(img);
dtype = class(img);
switch dtype
    case 'double'
        dtype = 'single';
        img = single(img);
    case 'logical'
        dtype = 'int8';
        img = int8(img);
end

% Adjust for saving to Nifti
orient = diag([-1 -1 1 1]) * orient;
d = d([2,1,3]);
fov = fov([2,1,3]);
img = permute(img,[2,1,3]);

info = init_niftiinfo(label,fov./d,dtype,d);
if gz_flag
    svname = fname(1:end-3);
else
    svname = fname;
end

info.TransformName = 'Sform';
info.Transform = affine3d(orient');
niftiwrite(img,svname,info,'Compressed',gz_flag);
status = exist(fname,'file');
