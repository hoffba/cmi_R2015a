function [img,label,fov,orient,info] = readNIFTI(fname,d)
%% reads Nifti *.nii file
[~, Nfile, ~] = fileparts(fname);

info = niftiinfo(fname);
% Check that dimensions match
if nargin==2 && ~isempty(d) && ~all(d(1:3)==info.ImageSize)
    warning('New dimensions (%s) do not match existing (%s)\n',...
        sprintf('%u %u %u',info.ImageSize),sprintf('%u %u %u',d(1:3)));
    img = []; label = {}; fov = []; orient = []; info = []; return;
end
img = double(niftiread(info));

fov = info.ImageSize.*info.PixelDimensions;

label = {Nfile};
d4 = size(img,4);
if d4>1
    label = strcat(label,'_',cellfun(@num2str,num2cell(1:d4),'UniformOutput',false));
end

orient = (info.Transform.T * diag([-1 -1 1 1]))';
    
