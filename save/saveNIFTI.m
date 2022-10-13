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
info = init_niftiinfo(label,fov./d(1:3),dtype,d);
gz_flag = strcmp(fname(end-2:end),'.gz');
if gz_flag
    svname = fname(1:end-3);
else
    svname = fname;
end
if ~isempty(orient) && ~all(all(orient==eye(4)))
    info.TransformName = 'Sform';
    info.Transform = affine3d((diag([-1 -1 1 1])*orient)');
%     tform = affine3d(orient);
%     if tform.isRigid
%         % Use qform (quaternions)
%     else
%         % Use sform (affine)
%     end
end
niftiwrite(img,svname,info,'Compressed',gz_flag);
status = exist(fname,'file');
