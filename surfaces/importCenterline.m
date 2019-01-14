function lineData = importCenterline(fname)
% Import centerline coordinates from Mimics .txt file
% Assumes only single branch set was generated

if exist(fname,'file')
    fid = fopen(fname,'r');
    if fid
        str = fread(fid,inf,'*char')';
        fclose(fid);
        
        % Find first branch definition:
        str = regexp(str,':\r\n\r\n\r\n(.*)\r\n\r\n','tokens');
        str = strsplit(str{1}{1},'\r\n');
        str = cellfun(@(x)strsplit(strtrim(x))',str,'UniformOutput',false);
        str = [str{:}]';
        labels = str(1,:);
        vals = str2double(str(2:end,:));
        
        % Add points to output structure:
        for i = 1:size(str,2)
            lineData.(labels{i}) = vals(:,i);
        end
    else
        error('Could not open file: %s',fname);
    end
end