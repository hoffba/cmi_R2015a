function p = readTransformParam(fname)

p=[];
if exist(fname,'file')
    if endsWith(fname,'.txt')
        
        % Read file:
        fid = fopen(fname,'rt');
        str = fread(fid,'*char')';
        fclose(fid);

        % Parse parameters into structure
        tok = regexp(str,'\(([A-Za-z]+) ([^\)]+)\)','tokens');
        for i = 1:numel(tok)
            tname = tok{i}{1};
            tval = tok{i}{2};
            tnum = str2num(tval);
            if ~isnan(tnum)
                % Convert to number
                tval = tnum;
            elseif strcmp(tval([1,end]),'""')
                % Strip quotes
                tval = tval(2:end-1);
            end
            p.(tname) = tval;
        end
        
    else
        warning('Input file needs to be a .txt'); return;
    end
else
    warning('Could not find file: %s',fname); return;
end