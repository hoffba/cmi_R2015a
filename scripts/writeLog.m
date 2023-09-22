function writeLog(fn,str,varargin)

% write to command window
fprintf(str,varargin{:});

% write to log file
if ischar(fn) && endsWith(fn,{'.txt','.log'})
    fid = fopen(fn,'a');
    if fid>0
        fprintf(fid,str,varargin{:});
        fclose(fid);
    else
        fprintf('Could not open log file.\n');
    end
end
