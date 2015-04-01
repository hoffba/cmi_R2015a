function stat = cmi_csvwrite(fname,A)
% Writes data in cell array A to CSV file

stat = false;
if isnumeric(A)
    A = num2cell(A);
end
if ischar(fname) && ~isempty(A) && iscell(A)
    tpath = fileparts(fname);
    if isdir(tpath)
        [nr,nc] = size(A);
        fid = fopen(fname,'w');
        if (fid ~= -1)
            for i = 1:nr
                for j = 1:nc
                    if ischar(A{i,j})
                        fmt = '%s';
                    else
                        fmt = '%d';
                    end
                    fprintf(fid,fmt,A{i,j});
                    if j==nc
                        if i~=nr
                            fprintf(fid,'\n');
                        end
                    else
                        fprintf(fid,',');
                    end
                end
            end
            fclose(fid);
            stat = true;
        end
    end
end