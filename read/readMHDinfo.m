function info = readMHDinfo(fname)

info = [];

fid = fopen(fname,'rb');
if fid>2

    % Read MHD header
    hstr = fread(fid,inf,'*char')';
    fclose(fid);

    % Parse image info
    tok = regexp(hstr,'(\w*) = ([^\n]*)','tokens');
    tok = vertcat(tok{:});

    % Convert numeric values
    for i = 1:size(tok,1)
        valstr = strsplit(tok{i,2});
        val = str2double(valstr);
        str_ind = isnan(val);
        [bool_ind,bool_b] = ismember(lower(valstr),{'true','false'});
        if all(bool_ind)
            val = bool_b==1;
        elseif any(str_ind)
            valstr{~str_ind} = val(~str_ind);
            valstr{bool_ind} = bool_b(bool_ind)==1;
            val = valstr;
        end
        tok{i,2} = val;
    end

    % Convert to info structure
    info = cell2struct(tok(:,2),tok(:,1));

else
    disp('Unable to read file: %s',fname);
end