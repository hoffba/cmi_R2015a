function [img,label,fov,orient,info] = readMHD(varargin)
% Reads .mhd and associated .raw files into the cmi program
img = [];

% Read info from .mhd file
fname = char(varargin{1});
if nargin==2
    origD = varargin{2};
else
    origD = [];
end

% Read MHD metadata
info = readMHDinfo(fname);
label = info.label;

[path,bname,~] = fileparts(fname);
if isfile(fullfile(path,[bname,'.raw']))
    rawfname = fullfile(path,[bname,'.raw']);
elseif isfile(fullfile(path,[bname,'.zraw']))
    rawfname = fullfile(path,[bname,'.zraw']);
elseif isfield(info.native_info,'ElementDataFile')
    rawfname = info.native_info.ElementDataFile;
end

% Image needs permutation to align with MiTAP geometry
perm = [2,1,3];
d_out = info.d(perm);

% Read in the .raw file
if ~isempty(origD) && ((length(origD)~=length(info.d)) || ~all(origD==d_out))
    warning('Dimensions do not match: current[%u %u %u] ~= new[%u %u %u]',origD,d_out);
elseif exist(rawfname,'file')
    fid = fopen(rawfname, 'r');
    if fid>2
        % Check if zipped:
        if isfield(info.native_info,'CompressedData') && info.native_info.CompressedData
            img = fread(fid,inf,'uchar=>uint8');
            import com.mathworks.mlwidgets.io.InterruptibleStreamCopier
            b = java.util.zip.InflaterInputStream(java.io.ByteArrayInputStream(img));
            isc = InterruptibleStreamCopier.getInterruptibleStreamCopier;
            c = java.io.ByteArrayOutputStream;
            isc.copyStream(b,c);
            img = double(typecast(c.toByteArray,info.Etype));
        else
            img = fread(fid,inf,info.Etype);
        end
        fclose(fid);
        if numel(img)~=prod([info.d,info.nv])
            img(prod([info.d,info.nv])) = 0; % in case file is incomplete we can see what's there
        end
        img = permute(reshape(img,[info.nv,info.d]),[2,3,4,1]);
    else
        disp('File could not be read correctly: %s',rawfname);
    end
end

img = permute(img,[perm,4]);
fov = info.fov(perm);
orient = info.orient([perm,4],[perm,4]);
info = info.native_info;

