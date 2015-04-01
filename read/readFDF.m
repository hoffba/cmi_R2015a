function [img,label,fov] = readFDF(varargin)
fnames = varargin{1};
[path,~] = fileparts(fnames);
% Find and load all FDF files in selected directory
fnames = dir([path filesep '*.fdf']);
% Read in each FDF in the directory
strarray = '';
valarray = '';
check = false;
for ifn = 1:size(fnames,1)
    % Read first FDF for header info
    fid = fopen([path filesep fnames(ifn).name],'r','l');
    num = 0;
    done = false;
    line = fgetl(fid);%disp(line)
    d4 = 1; ns = 1; echo = 1; vec = 1;
    while (~isempty(line) && ~done)
        line = fgetl(fid);%disp(line)
        % For first file, load general image details
        if ifn == 1
            if any(strncmp(line,{'float  rank = ','float rank = '},13))
                [rank, ~] = strtok(line,'float  rank = ;');
            end
            if any(strncmp(line,{'float  matrix[] = ','float matrix[] = '},17))
                [token, rem] = strtok(line,'float  matrix[] = { , };');
                d1 = str2num(token);
                d2 = str2num(strtok(rem,', };'));
                if strcmp(rank,'3')
                    rem_st = strfind(rem,',');
                    d3 = str2num(rem(rem_st(2)+1:end-2));
                else
                    d3 = 1;
                end
            end
            if any(strncmp(line,{'float  bits = ','float bits = '},13))
                [token, ~] = strtok(line,'float  bits = { , };');
                bits = str2num(token);
            end
            if strncmp('float  array_dim = ',line,19)
                [token, ~] = strtok(line,'float  array_dim = ');
                d4 = str2num(token);
            end
            if strncmp('int    slices = ',line,16)
                [token, ~] = strtok(line,'int    slices = ');
                ns = str2num(token);
            end
            if any(strncmp(line,{'float  roi[] = ','float roi[] = '},14))
                [token,rem] = strtok(line,'float  roi[] = ;');
                fov = cell2mat(textscan([token rem],'{%f,%f,%f}'));
            end
            if check % Read array label and current value
                if ~strncmp('int checksum = ',line,15)
                    token = textscan(line,'%s %s %s %f');
                    strarray = token{2}{1};
                    valarray = token{4};
                end
                check = false;
            end
            if strncmp('float  orientation[] = ',line,23)
                check = true; % Next line will be the array label
            end
        else
            if strncmp(['float  ' strarray ' = '],line,length(strarray)+10)
                valarray = str2num(strtok(line,['float  ' strarray ' = ;']));
            end
        end
        % After first image, only look for matrix index values
        if strncmp('int    slice_no = ',line,18)
            [token, ~] = strtok(line,'int    slice_no = ');
            slc = str2num(token);
        end
        if strncmp('int    echo_no = ',line,17)
            [token, ~] = strtok(line,'int    echo_no = ');
            echo = str2num(token);
        end
        if strncmp('int    array_index = ',line,21)
            [token, ~] = strtok(line,'int    array_index = ');
            vec = str2num(token);
        end
        
        num = num + 1;
        if num > 48
            done = true;
        end
    end
    % Initialize image matrix
    if ifn == 1
        img = zeros(d2,d1,d3*ns,d4);
        label = cell(1,d4);
        if strcmp(rank,'2')
            fov(3) = fov(3) * ns;
        end
        fov = fov*10;
        % Determine what echo to keep
        if strcmp(strarray,'gdiff') % sems_iadc, keep second echo (first is nav)
            techo = 2;
        else
            techo = 1;
        end
    end
    if echo == techo % For now, only load second echo (the image for sems_iadc)
        % Load FDF file
        fseek(fid, -d1*d2*d3*bits/8, 'eof');
                % If your image shifted, you may turn "shift...." off by 
                % add a "%" before fseek.
        img1 = fread(fid, 'float32');
        if d3 > 1
            slc = 1:d3;
        end
        img(:,:,slc,vec)=permute(reshape(img1,[d1 d2 d3]),[2 1 3]);
        if isempty(label{vec})
            label{vec} = [strarray '=' num2str(valarray)];
        end
    end
    fclose(fid);
end