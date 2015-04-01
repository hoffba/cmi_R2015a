function [img,label,fov] = readTIFF(varargin)
%% reads TIFF image data file (RGB or Stacked Index)

if nargin>0
    fname = varargin{1};
    if exist(fname,'file')
        info = imfinfo(fname);
        nh = info(1).Height;
        nw = info(1).Width;
        nf = length(info);
        ncolor = info(1).SamplesPerPixel;
        img = zeros(nh,nw,nf,ncolor);
        for i = 1:nf
            img(:,:,i,:) = permute(imread(fname,'TIFF','Index',i),[1,2,4,3]);
        end
        switch info(1).ColorType
            case 'grayscale'
                label = {'TIFF'};
            case 'truecolor'
                label = {'Red','Green','Blue'};
            otherwise
                label = repmat({'TIFF'},[1,ncolor]);
        end
        fov = [nw/info(1).XResolution,nh/info(1).YResolution,nf];
    end
end