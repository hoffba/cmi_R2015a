function stat = setParamValue(fname,fnout,varargin)

if isempty(fnout)
    fnout = fname;
end

stat = false;
if exist(fname,'file')
    if endsWith(fname,'.txt')
        nv = numel(varargin);
        if nv/2 == round(nv/2)
            
            % Read file:
            fid = fopen(fname,'rt');
            str = fread(fid);
            fclose(fid);
            
            % Loop over pairs:
            for i = 1:nv/2
                parname = varargin{2*i-1};
                parval = varargin{2*i};
                
                % Cast parameter value as a string:
                if isnumeric(parval)
                    parval = num2str(parval);
                elseif ischar(parval)
                    if ~startsWith(parval,'"')
                        parval = ['"',parval,'"'];
                    end
                else
                    warning('Input parameter value must be numeric or string.');
                end

                % Replace text
                if contains(str,['(',parname])
                    str = regexprep(str,['\(',parname,' [^\)]*\)'],...
                                        ['\(',parname,' ',parval,'\)']);
                else
                    str = [str,sprintf('\n(%s %s)\n',parname,parval)];
                end
            end
            
            % Write parameter file:
            fid = fopen(fnout,'wt');
            fwrite(fid,str);
            fclose(fid);

            stat = true;
        else
            warning('Inputs after file name must be Name/Value pairs'); return;
        end
    else
        warning('Input file needs to be a .txt'); return;
    end
else
    warning('Could not find file: %s',fname); return;
end