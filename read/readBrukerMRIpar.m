function p = readBrukerMRIpar(fname)

if ischar(fname)
    fname = {fname};
end
for i = 1:length(fname)
    [~,fstr] = fileparts(fname{i});
    if exist(fname{i},'file') && ismember(fstr,{'acqp','method','reco'})
        fid = fopen(fname{i},'r');
        while ~feof(fid)
            line = fgetl(fid);
            val = regexp(line,'##\$(\w*)=(.*)','tokens');
            if ~isempty(val)
                fld = val{1}{1};
                val = val{1}{2};
                sval = str2double(val);
                if ~isnan(sval)
                    % scalar value
                    val = sval;
                elseif strcmp(val(1),'(')
                    % array
                    if val(2)==' '
                        line = fgetl(fid);
                        if strcmp(line(1),'(')
                            % ignore array of arrays
                            val = [];
                        else
                            tval = sscanf(line,'%f');
                            if isempty(tval)
                                % string
                                val = line;
                            else
                                % numeric
                                n = sscanf(val(3:end-2),'%f, ')';
                                if length(n)==1
                                    n(2) = 1;
                                end
                                val = nan(n);
                                ct = length(tval);
                                val(1:ct) = tval;
                                while ct<prod(n)
                                    tval = sscanf(fgetl(fid),'%f');
                                    nval = length(tval);
                                    val(ct+(1:nval)) = tval;
                                    ct = ct+nval;
                                end
                                val = permute(val,[2,1,3:1:length(n)]);
                            end
                        end
                    else
                        % Ignore
                        val = [];
                    end
                end
                if ~isempty(val)
                    p.(fld) = val;
                end
            end
        end
        fclose(fid);
    else
        error('File not found: %s',fname{i})
    end
end