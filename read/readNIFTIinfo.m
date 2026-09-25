function info = readNIFTIinfo(fname)

[~, Nfile, ~] = fileparts(fname);
info = niftiinfo(fname);
info.format = 'nii';
info = struct('native_info',info,'label',char(Nfile));

% Adjust image info for cmi
info.d = info.native_info.ImageSize;
info.voxsz = info.native_info.PixelDimensions;
info.fov = info.d .* info.voxsz;
info.orient = diag([-1 -1 1 1]) * info.native_info.Transform.T';
