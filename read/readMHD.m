function [img,label,fov,info] = readMHD(varargin)
% Reads .mhd and associated .raw files into the cmi program
img = []; label = {}; fov = [];

% Read info from .mhd file
fname = varargin{1};
if nargin==2
    origD = varargin{2};
else
    origD = [];
end
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
    else
        Etype = 'Unknown';
    end
    if ~exist(rawfname,'file')
        ind = find(strcmp('ElementDataFile',hstr),1);
        if ~isempty(ind)
            rawfname = fullfile(path,hstr{ind+1});
        end
    end
    ind = find(strcmp('Labels',hstr),1);
    if isempty(ind)
        if nv>1
            label = strcat(bname,cellfun(@num2str,num2cell(1:nv)',...
                                         'UniformOutput',false))';
        else
            label = {bname};
        end
    else
        label = regexp(hstr{ind+1},'\"(.*?)\"','tokens');
        label = [label{:}]';
    end
    ind = find(strcmp('CompressedData',hstr),1);
    if isempty(ind)
        zp = false;
    else
        zp = strcmpi(hstr{ind+1},'true');
    end
%     ind = find(strcmp('ElementMin',hstr),1);
%     if ~isempty(ind)
%         emin = str2double(hstr{ind+1});
%     end
%     ind = find(strcmp('ElementMax',hstr),1);
%     if ~isempty(ind)
%         emax = str2double(hstr{ind+1});
%     end
    switch lower(Etype(5:end))
        case 'double'
            Etype = 'double';
        case 'float'
            Etype = 'single';
        case 'char'
            Etype = 'int8';
        case 'uchar'
            Etype = 'uint8';
        case 'short'
            Etype = 'int16';
        case 'ushort'
            Etype = 'uint16';
        case 'int'
            Etype = 'int32';
        case 'uint'
            Etype = 'uint32';
        otherwise
            Etype = 'single';
    end
    hchk = true;
else
    disp('File was not read, sonny');
end

% Read in the .raw file
if ~isempty(origD) && ((length(origD)~=length(d)) || ~all(origD==d([2,1,3])))
    warning('Dimensions do not match: current[%u %u %u] ~= new[%u %u %u]',origD,d([2,1,3]));
elseif hchk && exist(rawfname,'file')
    fid = fopen(rawfname, 'r');
    if fid>2
        % Check if zipped:
        if zp
            img = fread(fid,inf,'uchar=>uint8');
            import com.mathworks.mlwidgets.io.InterruptibleStreamCopier
            b = java.util.zip.InflaterInputStream(java.io.ByteArrayInputStream(img));
            isc = InterruptibleStreamCopier.getInterruptibleStreamCopier;
            c = java.io.ByteArrayOutputStream;
            isc.copyStream(b,c);
            img = double(typecast(c.toByteArray,Etype));
        else
            img = fread(fid,inf,Etype);
        end
        if numel(img)~=prod([d,nv])
            img(prod([d,nv])) = 0; % in case file is incomplete we can see what's there
        end
        img = permute(reshape(img,[nv,d]),[3,2,4,1]);
        fov = d.*voxsz;
        fclose(fid);
    else
        disp('File could not be read correctly. Heartbreaker.');
    end
end


