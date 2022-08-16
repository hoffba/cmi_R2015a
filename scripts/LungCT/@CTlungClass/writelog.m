function writelog(self,str,varargin)

% Write to command line:
fprintf(str,varargin{:});

% Write to log file:
fid = fopen(self.fn_log,'a');
if fid
    fprintf(fid,str,varargin{:});
    fclose(fid);
else
    fprintf('Could not open log file.\n');
end