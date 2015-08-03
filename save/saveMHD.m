function status = saveMHD(fname,img,labels,max_ext,~)
% Save as MHD file for Elastix
% * Automatically adds label to name if dim(4) > 1

% Permute from YXZ to XYZ:
img = permute(img,[2,1,3,4]);
max_ext = max_ext([2,1,3]);

[d(1),d(2),d(3),nv] = size(img);
[pathstr, fname, ~] = fileparts(fname);
if (nv>1)
    % Add label to file name
    if length(labels)<nv
        labels = cellfun(@num2str,num2cell(1:nv)');
    end
    ofnames = strcat(fname,'_',labels);
else
    ofnames = fname;
end
% Data types
if islogical(img)
    Etype = 'uint8';
else
    Etype = 'float';
end
element_types = struct('double','MET_DOUBLE',...
                       'int8',  'MET_CHAR',...
                       'uint8', 'MET_UCHAR',...
                       'float', 'MET_FLOAT',...
                       'int16', 'MET_SHORT',...
                       'uint16','MET_USHORT',...
                       'int32', 'MET_INT',...
                       'uint32','MET_UINT');
if ~isa(img,'double')
    img = double(img);
end
NDims = length(d);
% ElementNumberOfChannels = 1;
voxsz = max_ext./d;
if ischar(ofnames)
    ofnames = {ofnames};
end

for ifn = 1:nv
    % First check that the file name is correct
    mhfname = fullfile(pathstr,[ofnames{ifn},'.mhd']);
    rawfname = [ofnames{ifn},'.raw'];

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
    fprintf(fid,'ElementType = %s\n',element_types.(Etype));
    fprintf(fid,'ElementDataFile = %s\n',rawfname);
    fprintf(fid,'ElementNumberOfChannels = 1');
    % fprintf(fid,'ElementMin = %.6f',emin);
    % fprintf(fid,'ElementMax = %.6f',emax);
    fclose(fid);

    % write *.raw file
    fid = fopen(fullfile(pathstr,rawfname),'w');
    status = fwrite(fid,img(:,:,:,ifn),Etype);
    fclose(fid);
end

