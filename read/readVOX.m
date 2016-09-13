function [img,label,fov] = readVOX(varargin)

ffname = varargin{1};
origD = varargin{2};
[fdir,fname,ext] = fileparts(ffname);

% % Get image info:
% fid = fopen(fullfile(fdir,'DataInfo.xml'),'r');
% while ~feof(fid)
%     str = fgetl(fid);
%     
% end
% fclose(fid);

% Read image data:
inames = {'VolumeSize','VoxelSize','VolumeScale','Endian','Field'};
fid = fopen(ffname,'r');
go = true;
while go
    str = fgetl(fid);
    disp(str);
    [pname,rem] = strtok(str);
    rem = strtrim(rem);
    if ~isempty(rem)
        if any(strcmp(pname,inames))
            val = str2num(rem);
            if isempty(val)
                % keep as string
                val = rem;
            end
            info.(pname) = val;
        end
        if (length(rem)>1) && strcmp(rem(end-1:end),'##')
            go = false;
        end
    end
end
switch info.VoxelSize
    case 16
        fmt = 'int16';
    case 8
        fmt = 'int8';
end
switch info.Endian
    case 'L'
        mchn = 'l';
    case 'B'
        mchn = 'b';
end
img = fread(fid,prod(info.VolumeSize),fmt,0,mchn);
fclose(fid);

scl = regexp(info.Field,'Scale ([\d.-]*) Offset ([\d.-]*)','tokens');
scl = cellfun(@str2double,scl{1});
img = permute(reshape(img,info.VolumeSize),[2,1,3]) * scl(1) + scl(2);
fov = info.VolumeSize.*info.VolumeScale;
label = {fname};