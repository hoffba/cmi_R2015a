function varargout = getData(self,varargin)
% Retrieve data for use in processing

N = numel(varargin);
varargout = cell(N,1);
for i = 1:N
    dat_str = varargin{i};
    if isempty(self.dat.(dat_str).mat)
        fname = self.getFileName(dat_str);
        if exist(fname,'file')
% Load data from file
            writeLog(self.fn_log,'Reading %s from file: %s\n',dat_str,fname);
            info = niftiinfo(fname);
            if self.check_info(info,dat_str)
                self.dat.(dat_str).info = info;
                self.dat.(dat_str).mat = niftiread(info);
            else
                writeLog(self.fn_log,'Warning: info for %s does not match.\n',dat_str);
            end
        else
% Perform processing to generate requested data
            self.run(dat_str)
        end
    end
% Return loaded matrix
    varargout{i} = self.dat.(dat_str).mat;
end
