function info = readNIFTIhdr(fid)
% Reads .nii header
fseek(fid,0,'bof');
    %  Original header structures	
	%  struct header_key                     /* header key      */ 
	%       {                                /* off + size      */
	%       int sizeof_hdr                   /*  0 +  4         */
	%       char data_type[10];              /*  4 + 10         */
	%       char db_name[18];                /* 14 + 18         */
	%       int extents;                     /* 32 +  4         */
	%       short int session_error;         /* 36 +  2         */
	%       char regular;                    /* 38 +  1         */
	%       char hkey_un0;                   /* 39 +  1 */
	%       };                               /* total=40 bytes  */
	%
	% int sizeof_header   Should be 348.
	% char regular        Must be 'r' to indicate that all images and 
	%                     volumes are the same size. 
directchar = 'uchar=>char';
info.key.sizeof_hdr    = fread(fid,1,'int32')';
info.key.data_type     = deblank(fread(fid,10,directchar)');
info.key.db_name       = deblank(fread(fid,18,directchar)');
info.key.extents       = fread(fid, 1,'int32')';
info.key.session_error = fread(fid, 1,'int16')';
info.key.regular       = fread(fid, 1,directchar)';
info.key.hkey_un0      = fread(fid, 1,directchar)';
    %  Original header structures    
	%  struct image_dimension
	%       {                                /* off + size      */
	%       short int dim[8];                /* 0 + 16          */
    %       /*
    %           dim[0]      Number of dimensions in database; usually 4. 
    %           dim[1]      Image X dimension;  number of *pixels* in an image row. 
    %           dim[2]      Image Y dimension;  number of *pixel rows* in slice. 
    %           dim[3]      Volume Z dimension; number of *slices* in a volume. 
    %           dim[4]      Time points; number of volumes in database
    %       */
	%       float intent_p1;   % char vox_units[4];   /* 16 + 4       */
	%       float intent_p2;   % char cal_units[8];   /* 20 + 4       */
	%       float intent_p3;   % char cal_units[8];   /* 24 + 4       */
	%       short int intent_code;   % short int unused1;   /* 28 + 2 */
	%       short int datatype;              /* 30 + 2          */
	%       short int bitpix;                /* 32 + 2          */
	%       short int slice_start;   % short int dim_un0;   /* 34 + 2 */
	%       float pixdim[8];                 /* 36 + 32         */
	%	/*
	%		pixdim[] specifies the voxel dimensions:
	%		pixdim[1] - voxel width, mm
	%		pixdim[2] - voxel height, mm
	%		pixdim[3] - slice thickness, mm
	%		pixdim[4] - volume timing, in msec
	%					..etc
	%	*/
	%       float vox_offset;                /* 68 + 4          */
	%       float scl_slope;   % float roi_scale;     /* 72 + 4 */
	%       float scl_inter;   % float funused1;      /* 76 + 4 */
	%       short slice_end;   % float funused2;      /* 80 + 2 */
	%       char slice_code;   % float funused2;      /* 82 + 1 */
	%       char xyzt_units;   % float funused2;      /* 83 + 1 */
	%       float cal_max;                   /* 84 + 4          */
	%       float cal_min;                   /* 88 + 4          */
	%       float slice_duration;   % int compressed; /* 92 + 4 */
	%       float toffset;   % int verified;          /* 96 + 4 */
	%       int glmax;                       /* 100 + 4         */
	%       int glmin;                       /* 104 + 4         */
	%       };                               /* total=108 bytes */
info.dim.dim        = fread(fid,8,'int16')';
info.dim.intent_p1  = fread(fid,1,'float32')';
info.dim.intent_p2  = fread(fid,1,'float32')';
info.dim.intent_p3  = fread(fid,1,'float32')';
info.dim.intent_code = fread(fid,1,'int16')';
info.dim.datatype   = fread(fid,1,'int16')';
info.dim.bitpix     = fread(fid,1,'int16')';
info.dim.slice_start = fread(fid,1,'int16')';
info.dim.pixdim     = fread(fid,8,'float32')';
info.dim.vox_offset = fread(fid,1,'float32')';
info.dim.scl_slope  = fread(fid,1,'float32')';
info.dim.scl_inter  = fread(fid,1,'float32')';
info.dim.slice_end  = fread(fid,1,'int16')';
info.dim.slice_code = fread(fid,1,'uchar')';
info.dim.xyzt_units = fread(fid,1,'uchar')';
info.dim.cal_max    = fread(fid,1,'float32')';
info.dim.cal_min    = fread(fid,1,'float32')';
info.dim.slice_duration = fread(fid,1,'float32')';
info.dim.toffset    = fread(fid,1,'float32')';
info.dim.glmax      = fread(fid,1,'int32')';
info.dim.glmin      = fread(fid,1,'int32')';
%  Original header structures
	%  struct data_history       
	%       {                                /* off + size      */
	%       char descrip[80];                /* 0 + 80          */
	%       char aux_file[24];               /* 80 + 24         */
	%       short int qform_code;            /* 104 + 2         */
	%       short int sform_code;            /* 106 + 2         */
	%       float quatern_b;                 /* 108 + 4         */
	%       float quatern_c;                 /* 112 + 4         */
	%       float quatern_d;                 /* 116 + 4         */
	%       float qoffset_x;                 /* 120 + 4         */
	%       float qoffset_y;                 /* 124 + 4         */
	%       float qoffset_z;                 /* 128 + 4         */
	%       float srow_x[4];                 /* 132 + 16        */
	%       float srow_y[4];                 /* 148 + 16        */
	%       float srow_z[4];                 /* 164 + 16        */
	%       char intent_name[16];            /* 180 + 16        */
	%       char magic[4];   % int smin;     /* 196 + 4         */
	%       };                               /* total=200 bytes */
info.hist.descrip     = deblank(fread(fid,80,directchar)');
info.hist.aux_file    = deblank(fread(fid,24,directchar)');
info.hist.qform_code  = fread(fid,1,'int16')';
info.hist.sform_code  = fread(fid,1,'int16')';
info.hist.quatern_b   = fread(fid,1,'float32')';
info.hist.quatern_c   = fread(fid,1,'float32')';
info.hist.quatern_d   = fread(fid,1,'float32')';
info.hist.qoffset_x   = fread(fid,1,'float32')';
info.hist.qoffset_y   = fread(fid,1,'float32')';
info.hist.qoffset_z   = fread(fid,1,'float32')';
info.hist.srow_x      = fread(fid,4,'float32')';
info.hist.srow_y      = fread(fid,4,'float32')';
info.hist.srow_z      = fread(fid,4,'float32')';
info.hist.intent_name = deblank(fread(fid,16,directchar)');
info.hist.magic       = deblank(fread(fid,4,directchar)');

fseek(fid,253,'bof');
info.hist.originator  = fread(fid, 5,'int16')';

