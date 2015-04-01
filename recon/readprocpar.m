function pp = readprocpar(tdir,varargin)
% function for reading parameters from procpar file
% [par1,par2,...] = read_procpar('par1','par2',...);
% (input and output parameters must be in the same order)

if (nargin>1) && iscellstr(varargin)
    %varargout = cell(1,size(varargin,2));
    pp = struct;
    procfid = fopen(fullfile(tdir,'procpar'),'r','ieee-be');
    i = 1;
    while i <= size(varargin,2)
        line = fgetl(procfid);
        if ~ischar(line) % invalid file
            i = size(varargin,2)+1;
        else
            token = strtok(line);
            % Find matching variables
            k = find(strcmp(token,varargin),1);
            if ~isempty(k)
                % Next line contains values
                line = fgetl(procfid);
                [n,rem] = strtok(line);
                n = str2num(n);
                if rem(2) == '"'
                    pp.(varargin{k}) = rem(3:end-1);
                    if n > 1
                        for j = 1:n-1
                            lne = fgetl(procfid);
                            pp.(varargin{k}) = char(pp.(varargin{k}),lne(2:end-1));
                        end
                    end
                else
                    pp.(varargin{k}) = str2num(rem);
                end
                i=i+1;
            end
        end
    end
    fclose(procfid);
    ind = find(~isfield(pp,varargin));
    for i = 1:length(ind)
        pp.(varargin{ind(i)}) = [];
    end
else
    pp = [];
end


