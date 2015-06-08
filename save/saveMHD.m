function status = saveMHD(fname,img,~,max_ext,info)
% Save as MHD file for Elastix

if nargin<5
    info = [];
end
% First check that the file name is correct
[pathstr, fname, ~] = fileparts(fname);
mhfname = fullfile(pathstr,[fname,'.mhd']);
rawfname = [fname,'.raw'];

% Data types
if ~isa(img,'double')
    img = double(img);
end
imax = max(img(:));
imin = min(img(:));
flchk = ((std(img(:))/(imax-imin))<0.01);
if islogical(img)
    Etype = 'uint8';
elseif flchk
    Etype = 'float';
else
    if imin<0 % Signed
        imax = max(imax,-imin);
        if imax>2^31
            Etype = 'float';
        elseif imax>2^15
            Etype = 'int32';
        elseif imax>2^7
            Etype = 'int16';
        else
            Etype = 'int8';
        end
    else % Unsigned
        if imax>2^32
            Etype = 'float';
        elseif imax>2^16
            Etype = 'uint32';
        elseif imax>2^8
            Etype = 'uint16';
        else
            Etype = 'uint8';
        end
    end
end
Etype = 'float';
element_types = struct('double','MET_DOUBLE',...
                       'int8',  'MET_CHAR',...
                       'uint8', 'MET_UCHAR',...
                       'float', 'MET_FLOAT',...
                       'int16', 'MET_SHORT',...
                       'uint16','MET_USHORT',...
                       'int32', 'MET_INT',...
                       'uint32','MET_UINT');

d = ones(1,length(max_ext));
td = size(img);
d(1:length(td)) = td;
% Pos = max_ext./d;
NDims = length(d);
% ElementNumberOfChannels = 1;
voxsz = max_ext./d;

% write *.mhd file
fid = fopen(mhfname,'w');
fprintf(fid,'ObjectType = Image\n');
fprintf(fid,'NDims = %u\n',NDims);
fprintf(fid,'BinaryData = %s\n','True');
fprintf(fid,'BinaryDataByteOrderMSB = %s\n','False');
fprintf(fid,'CompressedData = %s\n','False');
fprintf(fid,'TransformMatrix = %u %u %u %u %u %u %u %u %u\n',[ 1 0 0 0 1 0 0 0 1 ]);
fprintf(fid,['Position =',repmat(' %.10f',1,NDims),'\n'],(1-d).*voxsz/2);
% fprintf(fid,['CenterOfRotation =',repmat(' %.10f',1,NDims),'\n'],zeros(1,NDims));
% fprintf(fid,['Position =',repmat(' %.10f',1,NDims),'\n'],voxsz/2);
% fprintf(fid,['CenterOfRotationPoint =',repmat(' %.10f',1,NDims),'\n'],d.*voxsz/2);
fprintf(fid,'AnatomicalOrientation = %s\n','RAI');
fprintf(fid,['ElementSpacing =',repmat(' %.6f',1,NDims),'\n'],voxsz);
fprintf(fid,['DimSize =',repmat(' %u',1,NDims),'\n'],d);
% fprintf(fid,'ElementNumberOfChannels = %u\n',ElementNumberOfChannels);
fprintf(fid,'ElementType = %s\n',element_types.(Etype));
fprintf(fid,'ElementDataFile = %s\n',rawfname);
% fprintf(fid,'ElementMin = %.6f',emin);
% fprintf(fid,'ElementMax = %.6f',emax);
fclose(fid);

% write *.raw file
fid = fopen(fullfile(pathstr,rawfname),'w');
status = fwrite(fid,img,Etype);
fclose(fid);


