function [kdata,info] = getVarianKdata(fname)

if nargin==0
    fname = uigetdir(pwd,'Select Varian Data Folder:');
end
if ischar(fname) && isdir(fname)
    if strcmp(fname(end-3:end),'.fid')
        fname = {fname};
    else
        fdir = fname;
        fname = dir([fdir,filesep,'*.fid']);
        fname = fullfile(fdir,{fname(:).name});
    end
elseif ~(iscellstr(fname) && all(cellfun(@(x)strcmp(x(end-3:end),'.fid'),fname)))
    error('Invalid input');
end

nf = length(fname);
for i = 1:nf
end