function [img,label,fov,orient,info] = readMHD(varargin)
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
if exist(fullfile(path,[bname,'.raw']),'file')
    rawfname = fullfile(path,[bname,'.raw']);
elseif exist(fullfile(path,[bname,'.zraw']),'file')
    rawfname = fullfile(path,[bname,'.zraw']);
end
hchk = false;
fid = fopen(fullfile(path,[bname,'.mhd']),'rb');
if fid>2
%     emin = nan; emax = nan;
    hstr = strtrim(strsplit(fread(fid,inf,'*char')',{'\n','='}));
    fclose(fid);
    
    hstr(cellfun(@isempty,hstr)) = [];
    istart = find(strcmp(hstr,'ObjectType'),1);
    info = struct(hstr{istart:end});
    flds = fieldnames(info);
    for i = 1:length(flds)
        orig_val = info.(flds{i});
        if strcmp(orig_val,'True')
            val = true;
        elseif strcmp(orig_val,'False')
            val = false;
        else
            val = str2num(orig_val);
            if isempty(val)
                val = orig_val;
            end
        end
        info.(flds{i}) = val;
    end

    d = info.DimSize;
    voxsz = info.ElementSpacing;
    if isfield(info,'ElementNumberOfChannels')
        nv = info.ElementNumberOfChannels;
    else
        nv = 1;
    end
    pos = zeros(1,3);
    if isfield(info,'Position')
        pos = info.Position;
    elseif isfield(info,'Offset')
        pos = info.Offset;
    end
    orient = [reshape(info.TransformMatrix,3,3)'*diag(voxsz),pos';0 0 0 1];
    Etype = info.ElementType;

    if ~exist(rawfname,'file')
        rawfname = info.ElementDataFile;
    end

    if isfield(info,'Labels')
        label = regexp(info.Labels,'\"(.*?)\"','tokens');
        label = [label{:}]';
    else
        if nv>1
            label = strcat(bname,cellfun(@num2str,num2cell(1:nv)',...
                                         'UniformOutput',false))';
        else
            label = {bname};
        end
    end

    if isfield(info,'CompressedData')
        zp = info.CompressedData;
    else
        zp = false;
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
        img = permute(reshape(img,[nv,d]),[2,3,4,1]);
        fov = d.*voxsz;
        fclose(fid);
    else
        disp('File could not be read correctly. Heartbreaker.');
    end
end

perm = [2,1,3];
img = permute(img,perm);
fov = fov(perm);
orient = orient([perm,4],[perm,4]);

