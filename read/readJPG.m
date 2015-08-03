function [img,label,fov] = readJPG(varargin)
img = []; label = {}; fov = [];

if nargin && ischar(varargin{1})
    fname = varargin{1};
else
    [fname,fpath] = uigetfile('*.bmp','Load BMP Images');
    fname = fullfile(fpath,fname);
end
if exist(fname,'file')
    [~,bname,~] = fileparts(fname);
    
    % Read image info and set FOV (pixels):
    info = imfinfo(fname);
    fov = [ info.Width , info.Height , 1 ];
    
    % Load image:
    img = double(permute(imread(fname),[1,2,4,3]));
    
    % Set label based on colormap:
    % ** assume RGB for now
    switch info.NumberOfSamples
        case 1
            label = {bname};
        case 3
            label = {'R','G','B'};
        otherwise
            label = cellfun(@num2str,num2cell((1:info.NumberOfSamples)'),...
                'UniformOutput',false);
    end
end