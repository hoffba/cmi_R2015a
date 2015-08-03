function status = saveMHD_multiCh(fname,img,labels,max_ext,~)
% Save as MHD file for Elastix

[d(1),d(2),d(3),nv] = size(img);
if length(labels)<nv
    labels = cellfun(@num2str,num2cell(1:5)');
end
[pathstr, fname, ~] = fileparts(fname);
% ofnames = strcat(fname,'_',labels);
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
d = d([2,1,3]);
max_ext = max_ext([2,1,3]);
voxsz = max_ext./d;

% First check that the file name is correct
mhfname = fullfile(pathstr,[fname,'.mhd']);
rawfname = [fname,'.raw'];

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

fprintf(fid,'ElementNumberOfChannels = %u\n',nv);
fprintf(fid,['Labels =',repmat(' "%s"',1,nv)],labels{:});
% fprintf(fid,'ElementMin = %.6f',emin);
% fprintf(fid,'ElementMax = %.6f',emax);
fclose(fid);

% write *.raw file
fid = fopen(fullfile(pathstr,rawfname),'w');
status = fwrite(fid,permute(img,[2,1,3,4]),Etype);
fclose(fid);

