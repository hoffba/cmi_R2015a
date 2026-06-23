function info = readMHDinfo(fname)

info = [];

fid = fopen(fname,'rb');
if fid>2

    % Read MHD header
    hstr = fread(fid,inf,'*char')';
    fclose(fid);

    % Parse image info
    tok = regexp(hstr,'(\w*) = ([^\n]*)','tokens');
    tok = strtrim(vertcat(tok{:}));

    % Convert numeric values
    for i = 1:size(tok,1)
        valstr = strsplit(tok{i,2});
        val = str2double(valstr);
        str_ind = isnan(val);
        [bool_ind,bool_b] = ismember(lower(valstr),{'true','false'});
        if all(bool_ind)
            val = bool_b==1;
        elseif any(str_ind)
            if ~all(str_ind)
                valstr{~str_ind} = val(~str_ind);
            end
            if any(bool_ind)
                valstr{bool_ind} = bool_b(bool_ind)==1;
            end
            val = valstr;
        end
        if iscell(val) && isscalar(val)
            val = val{1};
        end
        tok{i,2} = val;
    end

    % Convert to info structure
    info = cell2struct(tok(:,2),tok(:,1));

    % Parse native info into MiTAP format
    info = struct('native_info',info);
    info.native_info.format = 'mhd';
    info.d = info.native_info.DimSize;
    info.voxsz = ones(1,3);
    if isfield(info.native_info,'ElementSpacing')
        info.voxsz = info.native_info.ElementSpacing;
    end
    info.fov = info.d .* info.voxsz;
    if isfield(info,'ElementNumberOfChannels')
        info.nv = info.native_info.ElementNumberOfChannels;
    else
        info.nv = 1;
    end
    pos = zeros(1,3);
    if isfield(info.native_info,'Position')
        pos = info.native_info.Position;
    elseif isfield(info.native_info,'Offset')
        pos = info.native_info.Offset;
    end
    T = eye(3);
    if isfield(info.native_info,'TransformMatrix')
        T = reshape(info.native_info.TransformMatrix,3,3)';
    end
    info.orient = [T*diag(info.voxsz),pos';0 0 0 1];
    
    if isfield(info.native_info,'Labels')
        label = regexp(info.native_info.Labels,'\"(.*?)\"','tokens');
        info.label = [label{:}]';
    else
        [~,bname,~] = fileparts(char(fname));
        if info.nv>1
            info.label = strcat(bname,cellfun(@num2str,num2cell(1:info.nv)','UniformOutput',false))';
        else
            info.label = {bname};
        end
    end
       
    switch lower(info.native_info.ElementType(5:end))
        case 'double'
            info.Etype = 'double';
        case 'float'
            info.Etype = 'single';
        case 'char'
            info.Etype = 'int8';
        case 'uchar'
            info.Etype = 'uint8';
        case 'short'
            info.Etype = 'int16';
        case 'ushort'
            info.Etype = 'uint16';
        case 'int'
            info.Etype = 'int32';
        case 'uint'
            info.Etype = 'uint32';
    end

else
    disp('Unable to read file: %s',fname);
end