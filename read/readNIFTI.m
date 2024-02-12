function [img,label,fov,orient,info] = readNIFTI(fname,d)
%% reads Nifti *.nii file
[~, Nfile, ~] = fileparts(fname);

info = niftiinfo(fname);

% Adjust image info for cmi
info.format = 'nii';
dims = info.ImageSize([2,1,3]);
voxsz = info.PixelDimensions([2,1,3]);
fov = dims.*voxsz;
orient = diag([-1 -1 1 1]) * info.Transform.T';

% Check that dimensions match
if nargin==2 && ~isempty(d) && ~all(d(1:3)==dims)
    warning('New dimensions (%s) do not match existing (%s)\n',...
        sprintf('%u %u %u',dims),sprintf('%u %u %u',d(1:3)));
    img = []; label = {}; fov = []; orient = []; info = []; return;
end
img = double(niftiread(info));
img = permute(img,[2,1,3]);

label = {extractBefore(Nfile,'.nii')};
d4 = size(img,4);
if d4>1
    label = strcat(label,'_',cellfun(@num2str,num2cell(1:d4),'UniformOutput',false));
end
