function [img,label,fov,orient,info] = readNIFTI(varargin)
%% reads Nifti *.nii file
fname = varargin{1};
[~, Nfile, ~] = fileparts(fname);

info = niftiinfo(fname);
img = double(niftiread(info));

fov = info.ImageSize.*info.PixelDimensions;

label = {Nfile};
d4 = size(img,4);
if d4>1
    label = strcat(label,'_',cellfun(@num2str,num2cell(1:d4),'UniformOutput',false));
end

% Re-order for Matlab display
% img = permute(img,[2,1,3]);

orient = (info.Transform.T * diag([-1 -1 1 1]))';
    
