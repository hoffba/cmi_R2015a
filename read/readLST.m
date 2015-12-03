function t = readLST(fname)
% Read .lst file with timestamps for Bruker CT acquisition and PM
% Output: t = vector of times (seconds)

if nargin==0
    [fname,fpath] = uigetfile('*.lst','Select LST file:');
    if ischar(fname)
        fname = fullfile(fpath,fname);
    else
        return
    end
end

fid = fopen(fname,'r');
if fid
    str = strsplit(fread(fid,inf,'*char')',{'=','\r'});
    fclose(fid);
    t = cellfun(@(x)sscanf(x,'%02u:%02u:%02u:%03u')',str(5:2:end),...
        'UniformOutput',false);
    t = cellfun(@(x) 3600*x(1) + 60*x(2) + x(3) + x(4)/1000,t);
else
    error(['Could not open file: ' fname]);
end