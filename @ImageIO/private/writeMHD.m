function status = writeMHD(fname,img,info,oinfo)
% Save as MHD file for Elastix
% Syntax:
%   stat = writeMHD(fname,img,info)
%   stat = writeMHD(fname,img,info,oinfo)
% Inputs:
%   fname   = output file name
%   img     = image matrix (3D or 4D)
%   info    = CMI info
%   oinfo   = original image info
% * Automatically adds label to name if dim(4) > 1

[fpath,fname,~] = fileparts(fname);

% Set custom tag defaults:
defs = {'x_Institution',info.Institution;...
        'x_Manufacturer',info.Manufacturer;...
        'x_Model',info.Model;...
        'x_FieldStrength',info.FieldStrength;...
        'x_StationName',info.StationName;...
        'x_SerialNo',info.SerialNo;...
        'x_SoftwareVersion',info.SoftwareVer;...
        'x_StudyDescription',info.StudyDescription;...
        'x_SeriesDescription',info.SeriesDescription;...
        'x_ImageTypeString',info.ImageType;...
        'x_4thDimLabel',info.Label;...
        'x_ExternalComment_opt_',info.Comment};

% Permute from YX to XY:
img = permute(img,[2,1,3,4]);
if (nargin==4) && strcmp(oinfo.Format,'MHD')
    info = oinfo; % Use original info instead of CMI info
    flds = {'ElementSpacing','Offset','ElementSpacing',...
        'CenterOfRotation','DimSize',};
    for i = 1:length(flds)
        if isfield(info,flds{i})
            info.(flds{i})([1,2]) = info.(flds{i})([2,1]);
        end
    end
else
    info = struct('ObjectType',{'Image'},...
                  'NDims',{3},...
                  'BinaryData',{true},...
                  'BinaryDataByteOrderMSB',{false},...
                  'CompressedData',{false},...
                  'TransformMatrix',{info.Orientation},...
                  'Offset',{info.Origin([2,1,3])},...
                  'CenterOfRotation',{zeros(1,3)},...
                  'AnatomicalOrientation',{info.AnatomicalOrientation([2,1,3])},...
                  'ElementSpacing',{info.VoxelSpacing([2,1,3])},...
                  'ElementSize',{info.VoxelSize([2,1,3])},...
                  'DimSize',{info.Dim([2,1,3])},...
                  'ElementType',{''},...
                  'ElementDataFile',{''});
end

% Add TLC header tags:
for i = 1:size(defs,1)
    if ~isfield(info,defs{i,1})
        info.(defs{i,1}) = defs{i,2};
    end
end

nv = size(img,4);
if (nv>1)
    % Add label to file name
    if length(info.Labels)<nv
        labels = cellfun(@num2str,num2cell(1:nv)');
    else
        labels = info.Labels;
    end
    ofnames = strcat(fname,'_',labels);
else
    ofnames = {fname};
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
info.ElementType = element_types.(Etype);
if ~isa(img,'double')
    img = double(img);
end

legalchars = 'a-zA-Z0-9\-\_\.';
for ifn = 1:nv
    % First check that the file name is correct
    tname = regexprep(ofnames{ifn},['[^' legalchars ']'],'');
    mhdfname = fullfile(pathstr,[tname,'.mhd']);
    rawfname = [tname,'.raw'];

    % write *.mhd file
    flds = fieldnames(info);
    fid = fopen(mhdfname,'w');
    % Write TLC tags:
    fprintf(fid,'%% *****************************************************\n');
    fprintf(fid,'%% Institution = %s\n',info.x_Institution);
    fprintf(fid,'%% Manufacturer = %s\n',info.x_Manufacturer);
    fprintf(fid,'%% Model = %s\n',info.x_Model);
    fprintf(fid,'%% FieldStrength = %s\n',info.x_FieldStrength);
    fprintf(fid,'%% StationName = %s\n',info.x_StationName);
    fprintf(fid,'%% SerialNo = %s\n',info.x_SerialNo);
    fprintf(fid,'%% SoftwareVersion = %s\n',info.x_SoftwareVer);
    fprintf(fid,'%% StudyDescription = %s\n',info.x_StudyDescription);
    fprintf(fid,'%% SeriesDescription = %s\n',info.x_SeriesDescription);
    fprintf(fid,'%% ImageTypeString = %s\n',info.x_ImageType);
    fprintf(fid,'%% 4th dim label = %s\n',info.x_4thDimLabel);
    fprintf(fid,'%% ExternalComment(opt) = %s\n',info.x_ExternalComment_opt_);
    fprintf(fid,'%% *****************************************************');
    % Write normal MHD tags:
    for i = 1:find(strncmp('x_',flds,2))
        fprintf(fid,'%s = %s',flds{i},num2str(info.flds{i}));
    end
    fclose(fid);

    % write *.raw file
    fid = fopen(fullfile(fpath,rawfname),'w');
    status = fwrite(fid,img(:,:,:,ifn),Etype);
    fclose(fid);
end

