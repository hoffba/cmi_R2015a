function C = cmi_csvread(fname)
% Reads data from CSV file

C = {};
if ischar(fname) && exist(fname,'file')
    fid = fopen(fname,'r');
    if (fid >2)
        
        % Read all in as string
        str = fread(fid,'*char');
        fclose(fid);
        
        % Separate fields into cell array:
        str = strsplit(str','\n');
        nr = length(str);
        str = cellfun(@(x)strsplit(x,','),str,'UniformOutput',false);
        nc = cellfun(@length,str);
        C = cell(nr,max(nc));
        for i = 1:nr
            C(i,1:nc(i)) = str{i};
        end
        
        % Convert numeric fields into numbers:
        nC = cellfun(@str2double,C);
        nn = ~isnan(nC);
        C(nn) = nC(nn);
        
    end
end