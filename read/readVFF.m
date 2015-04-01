function [img,label,fov] = readVFF(varargin)
%% reads VFF data file (GEHC MicroView)
fid = fopen(varargin{1},'r','b');
num = 0;
done = false;
modality = [];
line = 1;
while (~isempty(line) && ~done)
    line = fgetl(fid);disp(line);
    if strncmp('rank=',line,5);
        rank = strtok(line,'rank=;');
    end
    if strncmp('size=',line,5)
        token = textscan(strtok(line,'size=;'),'%f %f %f');
        d1 = token{1};
        d2 = token{2};
        if strcmp(rank,'3')
            d3 = token{3};
        else
            d3 = 1;
        end
    end
    if strncmp('spacing=',line,8)
        token = textscan(strtok(line,'spacing=;'),'%f %f %f');
        if ~strcmp(rank,'3')
            token{3} = 1;
        end
        voxsz = cell2mat(token);
    end
    if strncmp('elementsize=',line,12)
        tsz = str2double(strtok(line,'elementsize=;'));
    end
    if strncmp('bits=',line,5)
        bits = str2num(strtok(line,'bits=;'));
    end
    if strncmp('modality=',line,9)
        modality = strtok(line,'modality=;');
    end
    if strncmp('format=',line,7)
        format = line(8:end-1);
    end
    if strncmp('gantryPosition=',line,16)
        label = {line(17:end-1)};
    end
    num = num + 1;
    if strcmp(line,char(12)) %num > 50
        done = true;
    end
end
if exist('voxsz','var') && exist('tsz','var')
    voxsz = voxsz * tsz;
end
fseek(fid, -d1*d2*d3*bits/8,'eof');
if isempty(modality) || strcmp(format,'base')
    img = fread(fid,['int' num2str(bits)]);
elseif strcmp(modality,'ST')
    img = fread(fid,['float' num2str(bits)]);
else
    img = fread(fid, ['int' num2str(bits)]);
end
% 
img = flipdim(permute(reshape(img,[d1,d2,d3]),[2,1,3]),1);
fov = [d1 d2 d3].*voxsz;
if ~exist('label','var')
    label = {'CT'};
end
fclose(fid);