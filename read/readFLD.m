function [img,label,fov] = readFLD(varargin)
%% reads FLD AVS pseudo-field file
fname = varargin{1};
parameters = {'ndim=','dim1=','dim2=','dim3=','nspace=','veclen=','data=','field=','label=','min_ext=','max_ext='};
numcheck = logical([1 1 1 1 1 1 0 0 0 1 1]); % check if parameter will be converted to number
pvals = cell(1,size(parameters,2)); % values of above-specified header variables to look for
maxlines = 30;                      % max # of lines to expect in the header
hend = false;                       % flags a found end of header
bytoff = 0;                         % for finding end of header 
count = 0;                          % current line read in
% labelLim = 100;                     % limit of # of char per label

fid = fopen(fname,'r','b');
if fid == -1 % Could not open the file
    error(['Could not upen the file: ' fname])
else
    disp(' ')
    disp(fname)
    % First, read in header info
    while (count < maxlines) && (~hend)
        line = fgetl(fid);
        if strncmp(line,char([12 12]),2)
            hend = true;
            bytoff = position + 2; % determines how far back to rewind before starting to read image in
        elseif ~strcmp(line(1),'#')
            disp(line);
            position = ftell(fid);
            ind = strncmp(line,parameters,5);
            ind = find(ind,1);
            if ~isempty(ind)
                eq = strfind(line,'=');
                pound = strfind(line,'#');
                if isempty(pound)
                    pound = length(line) + 1;
                end
                line = strtrim(line(eq(1)+1:pound(1)-1));
                if numcheck(ind) % needs conversion to number
                    pvals{ind} = str2num(line);
                elseif ind == find(strcmp('data=',parameters),1) % data type (short / float)
                    if strfind(line,'short')
                        pvals{ind} = 'short';
                    elseif strfind(line,'float')
                        pvals{ind} = 'float';
                    elseif strfind(line,'integer')
                        pvals{ind} = 'int16';
                    else
                        display('Warning, no valid data type was read. Assuming short.');
                        pvals{ind} = 'short';
                    end
                elseif ind == find(strcmp('label=',parameters),1) % labels
                    c = textscan(line,['%q']);
                    label = c{1}';
                else
                    pvals{ind} = strtrim(line); % trim blank space off
                end
            end
        end
        count = count + 1;
    end
    % if header end was found, read in image
    if hend
        fov = pvals{find(strcmp(parameters,'max_ext='),1)} ...
            - pvals{find(strcmp(parameters,'min_ext='),1)};
        d1 = pvals{find(strcmp(parameters,'dim1='),1)};
        d2 = pvals{find(strcmp(parameters,'dim2='),1)};
        d3 = pvals{find(strcmp(parameters,'dim3='),1)};
        d4 = pvals{find(strcmp(parameters,'veclen='),1)};
        dtype = pvals{find(strcmp(parameters,'data='),1)};
        if ~exist('label','var') || isempty(label) % label was not found in header
            label = cell(1,d4);
            label(:) = {fname};
        end
        if isempty(fov)
            fov = [d1 d2 d3];
        end
        label = label(1:d4);
    
        % Next, build image
        status = fseek(fid,bytoff,'bof');
        if status ~= -1
            img = fread(fid,d1*d2*d3*d4,dtype);
            img = reshape(img,[d4 d2 d1 d3]);
            img = permute(img,[3 2 4 1]);
        end
    else
        display('Could not find end of header! Aborting.');
    end
    fclose(fid);
end