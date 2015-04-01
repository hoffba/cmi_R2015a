function stat = saveANALYZE(fname,img,label,max_ext,~)
% Save as ANALYZE 7.5 file format
%  - .hdr file contains header info
%  - .img file contains image data

stat = false;

% Determine relevant image/file attributes:
[path,fnprefix] = fileparts(fname);
if islogical(img)
    precision = 'uchar';
    datatype = 2;
    bitpix = 8;
else
    precision = 'int16';
    datatype = 4;
    bitpix = 16;
end
[d(1),d(2),d(3),d(4)] = size(img);
voxsz = max_ext./d(1:3);
sclM = 1;
sclB = 0;

% Write .hdr file:
fid = fopen(fullfile(path,[fnprefix,'.hdr']),'w');
if fid>2
    
    % ****** Header Key:
    fwrite(fid,348,'int32'); % Size of header must be 348
    str = sprintf('% -10s','Analyze');
    fwrite(fid,str(1:10),'uchar'); % Data type
    str = sprintf('% -18s',date);
    fwrite(fid,str(1:18),'uchar'); % DB Name
    fwrite(fid,16384,'int32'); % Extents
    fwrite(fid,0,'int16'); % Session_Error
    fwrite(fid,'r','uchar'); % Regular
    fwrite(fid,'a','uchar'); % hkey_un0
    
    % ****** Image Dimensions
    % dim: [ #dimensions , x-dim , y-dim , z-dim , #volumes , undoc , undoc , undoc ]
    fwrite(fid,[4,d(2),d(1),d(3),d(4),0,0,0],'int16'); 
    str = sprintf('% -4s','mm');
    fwrite(fid,str(1:4),'uchar'); % vox_units
    str = sprintf('% -8s','');
    fwrite(fid,str(1:8),'uchar'); % cal_units
    fwrite(fid,0,'int16'); % Unused1
    fwrite(fid,datatype,'int16'); % DataType
    fwrite(fid,bitpix,'int16'); % bitpix
    fwrite(fid,0,'int16'); % dim_un0
    % pixdim: [ #dimensions , x-dim , y-dim , z-dim , #volumes , undoc , undoc , undoc ]
    fwrite(fid,[4,voxsz([2,1,3]),zeros(1,4)],'float'); % pixdim
    fwrite(fid,0,'float'); % vox_offset
    fwrite(fid,sclM,'float'); % Scale factor
    fwrite(fid,sclB,'float'); % zero-intercept
    fwrite(fid,0,'float'); % funused3
    fwrite(fid,0,'float'); % cal_max
    fwrite(fid,0,'float'); % cal_min
    fwrite(fid,0,'float'); % compressed
    fwrite(fid,0,'float'); % verified
    fwrite(fid,max(img(:)),'int32'); % glmax
    fwrite(fid,min(img(:)),'int32'); % glmin
    
    % ****** Data History (optional info)
    str = sprintf('"%s" ',label{:});
    str = sprintf('% -80s',str);
    fwrite(fid,str(1:80),'uchar'); % descrip
    str = sprintf('% -24s','');
    fwrite(fid,str(1:24),'uchar'); % aux_file
    fwrite(fid,0,'uchar'); % orient
    str = sprintf('% -10s','');
    fwrite(fid,str(1:10),'uchar'); % originator
    str = sprintf('% -10s','');
    fwrite(fid,str(1:10),'uchar'); % generated
    str = sprintf('% -10s','');
    fwrite(fid,str(1:10),'uchar'); % scannum
    str = sprintf('% -10s','');
    fwrite(fid,str(1:10),'uchar'); % patient_id
    str = sprintf('% -10s','');
    fwrite(fid,str(1:10),'uchar'); % exp_date
    str = sprintf('% -10s','');
    fwrite(fid,str(1:10),'uchar'); % exp_time
    str = sprintf('% -3s','');
    fwrite(fid,str(1:3),'uchar'); % hist_un0
    fwrite(fid,0,'int32'); % views
    fwrite(fid,0,'int32'); % vols_added
    fwrite(fid,0,'int32'); % start_field
    fwrite(fid,0,'int32'); % field_skip
    fwrite(fid,0,'int32'); % omax
    fwrite(fid,0,'int32'); % omin
    fwrite(fid,0,'int32'); % smax
    fwrite(fid,0,'int32'); % smin
    stat = fclose(fid);
end

% Write .img file:
if stat==0
    fid = fopen(fullfile(path,[fnprefix,'.img']),'w');
    if fid>2
        fwrite(fid,permute(img,[2,1,3,4]),precision);
        stat = fclose(fid);
    end
end





%     0 None                     (Unknown bit per voxel) % DT_NONE, DT_UNKNOWN 
%     1 Binary                         (ubit1, bitpix=1) % DT_BINARY 
%     2 Unsigned char         (uchar or uint8, bitpix=8) % DT_UINT8, NIFTI_TYPE_UINT8 
%     4 Signed short                  (int16, bitpix=16) % DT_INT16, NIFTI_TYPE_INT16 
%     8 Signed integer                (int32, bitpix=32) % DT_INT32, NIFTI_TYPE_INT32 
%    16 Floating point    (single or float32, bitpix=32) % DT_FLOAT32, NIFTI_TYPE_FLOAT32 
%    32 Complex, 2 float32      (Use float32, bitpix=64) % DT_COMPLEX64, NIFTI_TYPE_COMPLEX64
%    64 Double precision  (double or float64, bitpix=64) % DT_FLOAT64, NIFTI_TYPE_FLOAT64 
%   128 uint RGB                  (Use uint8, bitpix=24) % DT_RGB24, NIFTI_TYPE_RGB24 
%   256 Signed char            (schar or int8, bitpix=8) % DT_INT8, NIFTI_TYPE_INT8 
%   511 Single RGB              (Use float32, bitpix=96) % DT_RGB96, NIFTI_TYPE_RGB96
%   512 Unsigned short               (uint16, bitpix=16) % DT_UNINT16, NIFTI_TYPE_UNINT16 
%   768 Unsigned integer             (uint32, bitpix=32) % DT_UNINT32, NIFTI_TYPE_UNINT32 
%  1024 Signed long long              (int64, bitpix=64) % DT_INT64, NIFTI_TYPE_INT64
%  1280 Unsigned long long           (uint64, bitpix=64) % DT_UINT64, NIFTI_TYPE_UINT64 
%  1536 Long double, float128  (Unsupported, bitpix=128) % DT_FLOAT128, NIFTI_TYPE_FLOAT128 
%  1792 Complex128, 2 float64  (Use float64, bitpix=128) % DT_COMPLEX128, NIFTI_TYPE_COMPLEX128 
%  2048 Complex256, 2 float128 (Unsupported, bitpix=256) % DT_COMPLEX128, NIFTI_TYPE_COMPLEX128 
%
%    switch double(hdr.dime.datatype),
%    case   1,
%       hdr.dime.bitpix = int16(1 ); precision = 'ubit1';
%    case   2,
%       hdr.dime.bitpix = int16(8 ); precision = 'uint8';
%    case   4,
%       hdr.dime.bitpix = int16(16); precision = 'int16';
%    case   8,
%       hdr.dime.bitpix = int16(32); precision = 'int32';
%    case  16,
%       hdr.dime.bitpix = int16(32); precision = 'float32';
%    case  32,
%       hdr.dime.bitpix = int16(64); precision = 'float32';
%    case  64,
%       hdr.dime.bitpix = int16(64); precision = 'float64';
%    case 128,
%       hdr.dime.bitpix = int16(24); precision = 'uint8';
%    case 256 
%       hdr.dime.bitpix = int16(8 ); precision = 'int8';
%    case 511,
%       hdr.dime.bitpix = int16(96); precision = 'float32';
%    case 512 
%       hdr.dime.bitpix = int16(16); precision = 'uint16';
%    case 768 
%       hdr.dime.bitpix = int16(32); precision = 'uint32';
%    case 1024
%       hdr.dime.bitpix = int16(64); precision = 'int64';
%    case 1280
%       hdr.dime.bitpix = int16(64); precision = 'uint64';
%    case 1792,
%       hdr.dime.bitpix = int16(128); precision = 'float64';
%    otherwise
%       error('This datatype is not supported');
%    end




