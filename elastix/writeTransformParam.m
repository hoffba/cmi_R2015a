function stat = writeTransformParam(fname,p)

stat = false;

if ~endsWith(fname,'.txt')
    fname = strcat(fname,'.txt');
end
        
% Open file for writing
fid = fopen(fname,'wt');

if fid
    % Cast parameter structure into string:
    flds = fieldnames(p);
    nf = numel(flds);
    for i = 1:nf
        val = p.(flds{i});
        if isnumeric(val)
            val = num2str(val);
        elseif ischar(val)
            val = strcat('"',val,'"');
        end
        fprintf(fid,'(%s %s)\n',flds{i},val);
    end

    % Close file:
    fclose(fid);
    stat = true;
end
