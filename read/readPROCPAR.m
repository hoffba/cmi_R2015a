function pp = readPROCPAR(tdir,pstrs,fstr)
% Reads parameters from Varian procpar file
% Format: pp = readPROCPAR(tdir)
%         pp = readPROCPAR(tdir,cellstr)
% Inputs: tdir = directory containing procpar file
%         pstrs (optional) = cell array of strings of specific parameters
% Output: pp = structure with fields of all or specified parameters
%                   * specified pars not found are returned as empty fields

if nargin<3 || isempty(fstr) || ~ischar(fstr)
    fstr = 'procpar';
end

pp = struct;
procfid = fopen(fullfile(tdir,fstr),'r','ieee-be');
if procfid>0
    if (nargin<2) || ~iscellstr(pstrs)
        strchk = false;
    else
        strchk = true;
        for i = 1:length(pstrs)
            pp.(pstrs{i}) = [];
        end
    end
    go = true;
    while go
        line = fgetl(procfid);
        if ~ischar(line)
            go = false;
        else
            [str,rem] = strtok(line);
%             disp(str)
            dtype = str2num(['uint8(',strtok(rem),')']);
            line = fgetl(procfid);
            [nval,rem] = strtok(line);
            nval = str2double(nval);
            if ~strchk || any(strcmp(str,pstrs))
                switch dtype
                    case {2,4} % 'string' ; 'flag'
                        val = cell(nval,1);
                        val{1} = rem(3:end-1);
                        for i = 2:nval
                            line = fgetl(procfid);
                            val{i} = line(2:end-1);
                        end
                    otherwise
                        % {1,5} % 'real' ; 'frequency'
                        % {3,6} % 'delay' ; 'pulse'
                        % 7 % 'integer'
                        val = str2num(rem);
                end
                if ~isvarname(str)
%                     str = genvarname(str);
                    str = matlab.lang.makeValidName(str);
                end
                pp.(str) = val;
                fgetl(procfid);
            else
                if ~any(dtype==[2,4])
                    nval = 1;
                end
                for i = 1:nval
                    fgetl(procfid);
                end
            end
        end
    end
    fclose(procfid);
end
    
   

