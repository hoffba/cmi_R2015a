function [img,label,fov] = readMHD(varargin)
% Reads .mhd and associated .raw files into the cmi program
img = []; label = {}; fov = [];

% Read info from .mhd file
fname = varargin{1};
[path,bname,~] = fileparts(fname);
rawfname = fullfile(path,[bname,'.raw']);
hchk = false;
fid = fopen(fullfile(path,[bname,'.mhd']),'rb');
if fid>2
%     emin = nan; emax = nan;
    hstr = strtrim(strsplit(fread(fid,inf,'*char')',{'\n','='}));
    fclose(fid);
    ind = find(strcmp('DimSize',hstr),1);
    if ~isempty(ind)
        d = str2num(hstr{ind+1});
    end
    ind = find(strcmp('ElementSpacing',hstr),1);
    if ~isempty(ind)
        voxsz = str2num(hstr{ind+1});
    end
    ind = find(strcmp('ElementNumberOfChannels',hstr),1);
    if ~isempty(ind)
        nv = str2double(hstr{ind+1});
        if isnan(nv)
            nv = 1;
        end
    else
        nv = 1;
    end
    ind = find(strcmp('ElementType',hstr),1);
    if ~isempty(ind)
        Etype = hstr{ind+1};
    end
    if ~exist(rawfname,'file')
        ind = find(strcmp('ElementDataFile',hstr),1);
        if ~isempty(ind)
            rawfname = fullfile(path,hstr{ind+1});
        end
    end
%     ind = find(strcmp('ElementMin',hstr),1);
%     if ~isempty(ind)
%         emin = str2double(hstr{ind+1});
%     end
%     ind = find(strcmp('ElementMax',hstr),1);
%     if ~isempty(ind)
%         emax = str2double(hstr{ind+1});
%     end
    switch Etype
        case 'MET_DOUBLE'
            Etype = 'double';
        case 'MET_FLOAT'
            Etype = 'float';
        case 'MET_CHAR'
            Etype = 'int8';
        case 'MET_UCHAR'
            Etype = 'uint8';
        case 'MET_SHORT'
            Etype = 'int16';
        case 'MET_USHORT'
            Etype = 'uint16';
        case 'MET_INT'
            Etype = 'int32';
        case 'MET_UINT'
            Etype = 'uint32';
    end
    hchk = true;
else
    disp('File was not read, sonny');
end

% Read in the .raw file
if hchk && exist(rawfname,'file')
    fid = fopen(rawfname, 'r');
    if fid>2
        img = permute(reshape(fread(fid,inf,Etype),[nv,d]),[2,3,4,1]);
        fclose(fid);
        fov = d.*voxsz;
        if nv>1
            label = strcat(bname,cellfun(@num2str,num2cell(1:nv)',...
                                         'UniformOutput',false));
        else
            label = {bname};
        end
    else
        disp('File could not be read correctly. Heartbreaker.');
    end
end


