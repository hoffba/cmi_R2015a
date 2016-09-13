function status = saveMHD(fname,img,labels,fov,info)
% Save as MHD file for Elastix
% * Automatically adds label to name if dim(4) > 1

% Determine image info:
reqflds = {};
defval = {};
if nargin<5
    info = struct();
end
for i = 1:length(reqflds)
    if ~isfield(info,reqflds{i})
        info.(reqflds{i}) = defval{i};
    end
end

% Permute from YXZ to XYZ:
img = permute(img,[2,1,3,4]);
fov = fov([2,1,3]);

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
voxsz = fov./d;
if ischar(ofnames)
    ofnames = {ofnames};
end

legalchars = 'a-zA-Z0-9\-\_\.';
for ifn = 1:nv
    % First check that the file name is correct
    tname = regexprep(ofnames{ifn},['[^' legalchars ']'],'');
    mhfname = fullfile(pathstr,[tname,'.mhd']);
    rawfname = [tname,'.raw'];

    % write *.mhd file
    fid = fopen(mhfname,'w');
    
    % UMich header:
    fprintf(fid,'%% *****************************************************\n');
    
    
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

