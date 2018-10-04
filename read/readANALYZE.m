function [img,label,fov,info] = readANALYZE(varargin)
% Read ANALYZE 7.5 file format
%  - .hdr file contains header info
%  - .img file contains image data

img = []; label = {}; fov = [];
info = struct('SlicePos',{[0,0,0]},'SliceOrient',{[1,0,0,0,1,0]});

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
        info.hdrsz = fread(fid,1,'int32');
    end

    % ****** Header Key:
    info.hkDataType      = fread(fid,10,'uchar=>char');
    info.dbname          = fread(fid,18,'uchar=>char');
    info.hkExtent        = fread(fid,1 ,'int32'); % Should = 16384
    info.Session_Error   = fread(fid,1 ,'int16');
    info.hkreg           = fread(fid,1 ,'uchar=>char');
    info.hkey_un0        = fread(fid,1 ,'uchar=>char');

    % ****** Image Dimensions
    % dim: [ #dimensions , x-dim , y-dim , z-dim , #volumes , undoc , undoc , undoc ]
    info.dim             = fread(fid,8 ,'int16')';
    info.vox_units       = fread(fid,4 ,'uchar=>char')';
    info.cal_units       = fread(fid,8 ,'uchar=>char')';
    info.unused1         = fread(fid,1 ,'int16');
    info.DataType        = fread(fid,1 ,'int16');
    info.bitpix          = fread(fid,1 ,'int16');
    info.dim_un0         = fread(fid,1 ,'int16');
    info.voxsz           = fread(fid,8 ,'float')';
    info.vox_offset      = fread(fid,1 ,'float');
    info.sclM            = fread(fid,1 ,'float');
    info.sclB            = fread(fid,1 ,'float');
    info.funused3        = fread(fid,1 ,'float');
    info.cal_max         = fread(fid,1 ,'float');
    info.cal_min         = fread(fid,1 ,'float');
    info.compressed      = fread(fid,1 ,'float');
    info.verified        = fread(fid,1 ,'float');
    info.glmax           = fread(fid,1 ,'int32');
    info.glmin           = fread(fid,1 ,'int32');

    % ****** Data History (optional info)
    info.descrip         = fread(fid,80,'uchar=>char')';
    info.aux_file        = fread(fid,24,'uchar=>char')';
    info.orient          = fread(fid,1 ,'uchar');
    info.originator      = fread(fid,10,'uchar=>char')';
    info.generated       = fread(fid,10,'uchar=>char')';
    info.scannum         = fread(fid,10,'uchar=>char')';
    info.patient_id      = fread(fid,10,'uchar=>char')';
    info.exp_date        = fread(fid,10,'uchar=>char')';
    info.exp_time        = fread(fid,10,'uchar=>char')';
    info.hist_un0        = fread(fid,3 ,'uchar=>char')';
    info.views           = fread(fid,1,'int32');
    info.vols_added      = fread(fid,1,'int32');
    info.start_field     = fread(fid,1,'int32');
    info.field_skip      = fread(fid,1,'int32');
    info.omax            = fread(fid,1,'int32');
    info.omin            = fread(fid,1,'int32');
    info.smax            = fread(fid,1,'int32');
    info.smin            = fread(fid,1,'int32');
    stat = fclose(fid);
end

% Read .img file:
if stat==0
    estr = ''; 
    d = prod(info.dim(2:(info.dim(1)+1)));
    switch info.DataType
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
            estr = ['DataType = ',num2str(info.DataType)];
    end
    if isempty(estr)
        fid = fopen(imgname,'r',machstr);
        if fid>2
            img = fread(fid,d,precision);
            fclose(fid);
            if info.DataType==32 % complex
                img = complex(img(1:2:end),img(2:2:end));
            end
            img = reshape(img,info.dim(2:(info.dim(1)+1)));
            img = permute(img,[2,1,3,4]);
            fov = info.voxsz(2:4).*info.dim(2:4);
            label = strrep(regexp(info.descrip,'"(.[^"]*)"','match'),'"','');
            if length(label)~=info.dim(5)
                label = strcat('img',cellfun(@(x)sprintf('%02u',x),num2cell((1:info.dim(5))'),'UniformOutput',false)');
            end
        else
            error('IMG file not found!');
        end
    else
        error(['Invalid data type: ',estr]);
    end
end
