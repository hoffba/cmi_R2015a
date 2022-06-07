function writeLog(fn,str,varargin)

% write to command window
fprintf(str,varargin{:});

% write to log file
if ischar(fn) && exist(fn,'file') && endsWith(fn,'_log.txt')
    fid = fopen(fn,'a');
    if fid
        fprintf(fid,str,varargin{:});
        fclose(fid);
    else
        fprintf('Could not open log file.\n');
    end
end
