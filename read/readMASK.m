function [img,label,fov] = readMASK(varargin)
%% read in .mask type file as VOI
fname = varargin{1};
d = varargin{2}(1:3);
fid = fopen(fname);
if fid~=-1
    img = fread(fid,inf,'ubit1');
%     if length(img) < prod(d(1:3))
%         error('ABORT: Selected mask is smaller than the current image!')
%     elseif length(img) > prod(d(1:3))
%         error('ABORT: Selected mask is larger than the current image!')
%     else
%         img = reshape(img,d);
%     end
    nel = length(img); nout = prod(d);
    if nel < nout
        img((nel+1):nout) = 0;
    end
    img = reshape(img(1:nout),d);
else
    error('Could not open selected file. Check permissions.');
end
fclose(fid);
fov = [];
label = [];