function [img,label,fov] = readFDF(varargin)
fnames = varargin{1};
[path,~] = fileparts(fnames);
% Find and load all FDF files in selected directory
fnames = dir([path filesep '*.fdf']);
nf = length(fnames);
% Read in each FDF in the directory
strarray = '';
valarray = '';
for ifn = 1:nf
    % Read first FDF for header info
    fid = fopen([path filesep fnames(ifn).name],'r','l');
    done = false;
    line = fgetl(fid);%disp(line)
    na = 1; vec = 1; ns = 1; ne = 1; echo = 1;
    while (~isempty(line) && ~done)
        
        line = fgetl(fid);%disp(line)
        line = strsplit(line,' = ');
        
        if length(line)==2
            svar = strsplit(line{1});
            svar = svar{end};
            sval = strtok(line{end},';');
            
            if strcmp(svar(1),'*')
                % remove * from variable name
                svar(1) = [];
                % find strings between quotes
                sval = regexp(sval,'(?<=")[^"]+(?=")','match');
                sval = sval(1:2:end);
            elseif strcmp('[]',svar(end-1:end))
                % remove [] from variable name
                svar(end-1:end) = [];
                % find numbers in {} separated by comma
                sval = cellfun(@str2double,strsplit(sval(2:end-1),','));
            else
                % sval is a scalar
                sval = str2double(sval);
            end
            
        else
            svar = '';
            sval = '';
        end
        
        if ifn==1
        % For first file, load general image details
            switch svar
                case 'rank'
                    rank = sval;
                case 'matrix'
                    d1 = sval(1);
                    d2 = sval(2);
                    if rank==3
                        d3 = sval(3);
                    else
                        d3 = 1;
                    end
                case 'echoes'
                    ne = sval;
                case 'bits'
                    bits = sval;
                case 'array_dim'
                    na = sval;
                case 'slices'
                    ns = sval;
                case 'roi'
                    fov = sval;
                case 'array_name'
                    if ~strcmp(sval,'none')
                        strarray = svar;
                    end
            end
        end
        % After first image, only look for matrix index values
        switch svar
            case 'slice_no'
                slc = sval;
            case 'echo_no'
                echo = sval;
            case 'array_index'
                vec = sval;
            case strarray
                valarray = sval;
        end
        
        if strcmp(svar,'checksum')
            done = true;
        end
    end
    % Initialize image matrix
    if ifn == 1
        img = zeros(d2,d1,d3*ns,ne,na);
        label = cell(1,na);
        if rank==2
            fov(3) = fov(3) * ns;
        end
        fov = fov*10; % convert to mm
    end
    if ~(strcmp(strarray,'gdiff') && (echo==1))
        % *For sems_iadc, keep only second echo (first is nav)
        %   otherwise, load all echoes
        
        % Load FDF file
        fseek(fid, -d1*d2*d3*bits/8, 'eof');
                % If your image shifted, you may turn "shift...." off by 
                % add a "%" before fseek.
        img1 = fread(fid, 'float32');
        if d3 > 1
            slc = 1:d3;
        end
        img(:,:,slc,echo,vec)=permute(reshape(img1,[d1 d2 d3]),[2 1 3]);
        if isempty(label{vec})
            label{vec} = [strarray '=' num2str(valarray)];
        end
    end
    fclose(fid);
end
img = reshape(img,d1,d2,d3*ns,[]);
label = cellfun(@num2str,num2cell(1:size(img,4)),'UniformOutput',false);
