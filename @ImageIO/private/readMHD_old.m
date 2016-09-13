function [img,info] = readMHD(fname,D)
% Syntax:
%   [img,info] = readMHD(fname)
%   [img,info] = readMHD(fname,D)
%       - D specifies desired image dimensions [1x3];
%           load fails if dimensions don't match
% Reads .mhd and associated .raw files into the cmi program
img = [];
info = [];

% Read info from .mhd file
[path,bname,~] = fileparts(fname);
rawfname = fullfile(path,[bname,'.raw']);
fid = fopen(fullfile(path,[bname,'.mhd']),'rb');
if fid>2
    
    % Read image info from .mhd head file:
    info = struct('FileName',{fname},...
                  'Format',{'MHD'});
    lstr = fgetl(fid);
    while ischar(lstr)
        lstr = regexp(lstr,'(.*) = (.*)','tokens');
        if length(lstr)==1 %Only finds one '='
            fld = lstr{1}{1};
            val = lstr{1}{2};
        % Validate fieldname string:
            fld = matlab.lang.makeValidName(fld); 
            fld = matlab.lang.makeUniqueStrings(fld,fieldnames(info));
        % Convert numeric values
            nval = str2num(val);
            if ~isempty(nval)
                val = nval;
            end
        % Set info value
            info.(fld) = val;
        end
        lstr = fgetl(fid);
    end
    fclose(fid);
    
    % Check dimension size match
    nD = ones(1,3);
    nD(1:length(info.DimSize)) = info.DimSize;
    if (nargin==2) && ~isempty(D) && ~all(nD==D)
        return;
    end
    
    % Determine bitdepth
    switch info.ElementType
        case 'MET_DOUBLE'
            Etype = 'double';
        case 'MET_FLOAT'
            Etype = 'float';
        case 'MET_CHAR'
            Etype = 'int8';
        case 'MET_UCHAR'
            Etype = 'uint8';
        case 'MET_SHORT'
            Etype = 'int16';
        case 'MET_USHORT'
            Etype = 'uint16';
        case 'MET_INT'
            Etype = 'int32';
        case 'MET_UINT'
            Etype = 'uint32';
    end
    
    % 4D images:
    if isfield(info,'ElementNumberOfChannels')
        nv = info.ElementNumberOfChannels;
    else
        nv = 1;
    end
    nD = [nD,nv];
    
    % Read in the .raw file
    if exist(rawfname,'file')
        fid = fopen(rawfname, 'r');
        if fid>2
            img = fread(fid,inf,Etype);
            if length(img)<prod(nD)
                % in case file is incomplete we can see what's there
                img(prod(nD)) = 0; 
            end
            img = permute(reshape(img,nD),[2,1,3,4]);
            fclose(fid);
        else
            disp('File could not be opened. Check permissions.');
        end
    end
    
    % Permute geometry info from XY to YX:
    flds = {'ElementSpacing','Offset','ElementSpacing',...
        'CenterOfRotation','DimSize',};
    for i = 1:length(flds)
        if isfield(info,flds{i})
            info.(flds{i})([1,2]) = info.(flds{i})([2,1]);
        end
    end
end



