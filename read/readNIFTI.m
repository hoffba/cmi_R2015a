function [img,label,fov,orient,info] = readNIFTI(fname,d)
%% reads Nifti *.nii file
[~, Nfile, ~] = fileparts(fname);
label = {Nfile};

info = readNIFTIinfo(fname);

% Needs permutation to align with MiTAP geometry
perm = [2,1,3];
d_out = info.d(perm);
fov = info.fov(perm);
orient = info.orient;
info = info.native_info;

% Check that dimensions match
if nargin==2 && ~isempty(d) && ~all(d(1:3)==d_out)
    warning('New dimensions (%s) do not match existing (%s)\n',sprintf('%u %u %u',d_out),sprintf('%u %u %u',d(1:3)));
    img = []; label = {}; fov = []; orient = []; info = []; return;
end
img = double(niftiread(info));
img = permute(img,[perm,5,4]);

d4 = size(img,4);
if d4>1
    label = strcat(label{1},'_',cellfun(@num2str,num2cell(1:d4),'UniformOutput',false));
end
