function varargout = cmi_readPROCPAR(fdir,pars)
% function for reading parameters from procpar file
% [par1,par2,...] = read_procpar('par1','par2',...);
% (input and output parameters must be in the same order)

varargout = cell(1,length(pars));
procfid = fopen(fullfile(fdir,'procpar'),'r','ieee-be');
i = 1;
line = 'a';
% Read each line of the "procpar" file until either all parameters are 
% found or the search reaches the end of the file.
while (i <= length(pars) && ischar(line))
    line = fgetl(procfid);
    %disp(line);
    if ischar(line) % Make sure we are not at end of file (EOF)
        token = strtok(line); % Grab first word
        k = find(strcmp(token,pars)); % Find parameter index if requested
        if ~isempty(k)
            line = fgetl(procfid); % Parameter
            [n,rem] = strtok(line);
            n = str2num(n); % First is number of values for this parameter
            if rem(2) == '"' % If string
                varargout{k} = rem(3:end-1);
                if n > 1
                    for j = 1:n-1
                        % in loop to remove quotes
                        lne = fgetl(procfid);
                        varargout{k} = char(varargout{k},lne(2:end-1));
                    end
                end
            else % If number
                varargout{k} = str2num(rem);
            end
            i=i+1;
        end
    end
end
fclose(procfid);

