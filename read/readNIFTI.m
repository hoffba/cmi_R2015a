function [img,label,fov] = readNIFTI(varargin)
%% reads Nifti *.nii file
fname = varargin{1};
[~, Nfile, ~] = fileparts(fname);

info = niftiinfo(fname);
img = double(niftiread(info));

fov = info.ImageSize.*info.PixelDimensions;

if size(img,4)>1
    for i = 1:size(img,4)
        label(i) = {[Nfile,'_',num2str(i)]};
    end
else
    label = {Nfile};
end
    
