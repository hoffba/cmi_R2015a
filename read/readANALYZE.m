function [img,label,fov] = readANALYZE(varargin)
% Read ANALYZE 7.5 file format
%  - .hdr file contains header info
%  - .img file contains image data

img = []; label = {}; fov = [];

% Determine relevant image/file attributes:
fname = varargin{1};
[path,bname] = fileparts(fname);
hdrname = fullfile(path,[bname,'.hdr']);
imgname = fullfile(path,[bname,'.img']);

% Test endian:
machstr = 'l';
fid = fopen(hdrname,'r',machstr);
if fid>2
    % Read .hdr file:
    % Format information found at:
    % http://www.grahamwideman.com/gw/brain/analyze/formatdoc.htm
    
    % Check for correct endian:
    hdrsz = fread(fid,1,'int32');
    if hdrsz ~= 348
        fclose(fid);
        machstr = 'b';
        fid = fopen(hdrname,'r',machstr);
        hdrsz = fread(fid,1,'int32');
    end

    % ****** Header Key:
    hkDataType      = fread(fid,10,'uchar=>char');
    dbname          = fread(fid,18,'uchar=>char');
    hkExtent        = fread(fid,1 ,'int32'); % Should = 16384
    Session_Error   = fread(fid,1 ,'int16');
    hkreg           = fread(fid,1 ,'uchar=>char');
    hkey_un0        = fread(fid,1 ,'uchar=>char');

    % ****** Image Dimensions
    % dim: [ #dimensions , x-dim , y-dim , z-dim , #volumes , undoc , undoc , undoc ]
    dim             = fread(fid,8 ,'int16')';
    vox_units       = fread(fid,4 ,'uchar=>char')';
    cal_units       = fread(fid,8 ,'uchar=>char')';
    unused1         = fread(fid,1 ,'int16');
    DataType        = fread(fid,1 ,'int16');
    bitpix          = fread(fid,1 ,'int16');
    dim_un0         = fread(fid,1 ,'int16');
    voxsz           = fread(fid,8 ,'float')';
    vox_offset      = fread(fid,1 ,'float');
    sclM            = fread(fid,1 ,'float');
    sclB            = fread(fid,1 ,'float');
    funused3        = fread(fid,1 ,'float');
    cal_max         = fread(fid,1 ,'float');
    cal_min         = fread(fid,1 ,'float');
    compressed      = fread(fid,1 ,'float');
    verified        = fread(fid,1 ,'float');
    glmax           = fread(fid,1 ,'int32');
    glmin           = fread(fid,1 ,'int32');

    % ****** Data History (optional info)
    descrip         = fread(fid,80,'uchar=>char')';
    aux_file        = fread(fid,24,'uchar=>char')';
    orient          = fread(fid,1 ,'uchar');
    originator      = fread(fid,10,'uchar=>char')';
    generated       = fread(fid,10,'uchar=>char')';
    scannum         = fread(fid,10,'uchar=>char')';
    patient_id      = fread(fid,10,'uchar=>char')';
    exp_date        = fread(fid,10,'uchar=>char')';
    exp_time        = fread(fid,10,'uchar=>char')';
    hist_un0        = fread(fid,3 ,'uchar=>char')';
    views           = fread(fid,1,'int32');
    vols_added      = fread(fid,1,'int32');
    start_field     = fread(fid,1,'int32');
    field_skip      = fread(fid,1,'int32');
    omax            = fread(fid,1,'int32');
    omin            = fread(fid,1,'int32');
    smax            = fread(fid,1,'int32');
    smin            = fread(fid,1,'int32');
    stat = fclose(fid);
end

% Read .img file:
if stat==0
    estr = ''; 
    d = prod(dim(2:(dim(1)+1)));
    switch DataType
        case 0 % DT_NONE or DT_UNKNOWN
            estr = 'DT_NONE / DT_UNKNOWN';
        case 1 % DT_BINARY
            precision = 'ubit1';
        case 2 % DT_UNSIGNED_CHAR
            precision = 'uchar';
        case 4 % DT_SIGNED_SHORT
            precision = 'int16';
        case 8 % DT_SIGNED_INT
            precision = 'int32';
        case 16 % DT_FLOAT
            precision = 'float';
        case 32 % DT_COMPLEX
            precision = 'single';
            d = 2*d;
        case 64 % DT_DOUBLE
            precision = 'double';
        case 128 % DT_RGB
            estr = 'DT_RGB';
        case 255 % DT_ALL
            estr = 'DT_ALL';
        otherwise
            estr = ['DataType = ',num2str(DataType)];
    end
    if isempty(estr)
        fid = fopen(imgname,'r',machstr);
        if fid>2
            img = fread(fid,d,precision);
            fclose(fid);
            if DataType==32 % complex
                img = complex(img(1:2:end),img(2:2:end));
            end
            img = reshape(img,dim(2:(dim(1)+1)));
            img = permute(img,[2,1,3,4]);
            fov = voxsz(2:4).*dim(2:4);
            label = strrep(regexp(descrip,'"(.[^"]*)"','match'),'"','');
            if length(label)~=dim(5)
                label = strcat('img',cellfun(@(x)sprintf('%02u',x),num2cell((1:dim(5))'),'UniformOutput',false)');
            end
        else
            error('IMG file not found!');
        end
    else
        error(['Invalid data type: ',estr]);
    end
end
