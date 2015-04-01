% Function to read Elastix parameters from .txt file
function s = readElxTxt2Struct(fname)

if (nargin==0)
    [fname,path] = uigetfile('*.txt','Load ElastixParameters.txt',...
        'ElastixParameters.txt');
    if fname
        fname = fullfile(path,fname);
    else
        fname = '';
    end
elseif ~ischar(fname)
    error('Invalid input: fname');
end
if ~isempty(fname) && exist(fname,'file')
    % Read entire .txt file at once:
    fid = fopen(fname,'r');
    if fid>2
        str = fread(fid,'*char')';
        fclose(fid);
    else
        str = '';
    end
    % Parse variables into structure:
    if ~isempty(str)
        % Separate by line
        A = strtrim(strsplit(str,'\n'));
        for i = 1:length(A)
            % only look for statements in ()
            if ~isempty(A{i}) && strcmp(A{i}([1,end]),'()')
                [vname,val] = strtok(A{i}(2:end-1));
                if strcmp(val(2),'"') % "" indicate strings
                    q = strfind(val,'"');
                    nq = length(q);
                    if mod(nq,2) == 0
                        cstr = cell(1,nq/2);
                        for j = 1:nq/2
                            cstr{j} = val((q(2*j-1)+1):(q(2*j)-1));
                        end
                        if nq==2
                            cstr = cstr{1};
                        end
                        val = cstr;
                    else
                        warning(['Invalid string variable: ',vname])
                        val = [];
                    end
                else % otherwise assume numbers
                    val = str2num(val);
                end
                if ~isempty(val)
                    s.(vname) = val;
                end
            end
        end
    end
end

