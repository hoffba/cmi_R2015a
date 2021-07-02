function status = saveMASK(fname,img,varargin)
%% Save .mask file (binary VOI data)
% First check that the file name is correct
[pathstr, name, ext] = fileparts(fname);
if (isempty(ext) || ~strcmp(ext,'.mask'))
    fname = [pathstr name '.mask'];
end
fidmask = fopen(fname,'w');
count = fwrite(fidmask,img,'ubit1');
status = logical(count);