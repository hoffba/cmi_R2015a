function stat = file(self,method,dat_str)

N = numel(dat_str);
stat = false(1,N);

% Generate file names
fn_save = cellfun(@(x)sprintf('%s.%s.nii.gz',self.fn_base,self.dat.(x).tag),dat_str);

for i = 1:N
    switch method
        case 'check'
            stat(i) = logical(exist(fullfile(self.procdir,fn_save{i}),'file'));
        case 'save'
            if ~isempty(self.dat.(dat_str{i}).mat)
                try
                    niftiwrite(self.dat.(dat_str{i}).mat,fn_save{i}(1:end-3),self.dat.(dat_str{i}).info,'Compressed',true);
                    stat(i) = true;
                catch
                    writeLog(self.fn_log,'WARNING: Could not write file: %s\n',fn_save{i});
                end
            end
        case 'load'
            if exist(fullfile(self.procdir,fn_save{i}),'file')
                
            else
            end
        otherwise
    end
end