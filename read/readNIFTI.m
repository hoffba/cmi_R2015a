function [img,label,fov] = readNIFTI(varargin)
%% reads Nifti *.nii file
fname = varargin{1};

fid = fopen(fname,'r','ieee-le');
if fid>2
    info = readNIFTIhdr(fid);
    [info.location,info.fname] = fileparts(fname);
    switch info.dim.datatype
        case   1,
            info.dim.bitpix = 1;  precision = 'ubit1';
        case   2,
            info.dim.bitpix = 8;  precision = 'uint8';
        case   4,
            info.dim.bitpix = 16; precision = 'int16';
        case   8,
            info.dim.bitpix = 32; precision = 'int32';
        case  16,
            info.dim.bitpix = 32; precision = 'float32';
        case  32,
            info.dim.bitpix = 64; precision = 'float32';
        case  64,
            info.dim.bitpix = 64; precision = 'float64';
        case 128,
            info.dim.bitpix = 24; precision = 'uint8';
        case 256 
            info.dim.bitpix = 8;  precision = 'int8';
        case 511 
            info.dim.bitpix = 96; precision = 'float32';
        case 512 
            info.dim.bitpix = 16; precision = 'uint16';
        case 768 
            info.dim.bitpix = 32; precision = 'uint32';
        case 1024
            info.dim.bitpix = 64; precision = 'int64';
        case 1280
            info.dim.bitpix = 64; precision = 'uint64';
        case 1792,
            info.dim.bitpix = 128; precision = 'float64';
        otherwise
            error('This datatype is not supported'); 
    end
    fseek(fid,info.dim.vox_offset,'bof');
    img = reshape(fread(fid,inf,precision),info.dim.dim(2:4));
    fclose(fid);
    dorder = [3,2,4];
    fov = info.dim.dim(dorder).*info.dim.pixdim(dorder);
    label = {info.fname};
end



