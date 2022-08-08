function fname = getFileName(self,dat_str)

re_str = 're_'; % File name flag for resampled data
cell_flag = iscell(dat_str);
if ~cell_flag
    dat_str = {dat_str};
end
fname = cell(numel(dat_str),1);
for i = 1:numel(dat_str)
    if startsWith(dat_str{i},re_str)
        dat_str{i}(1:3) = [];
    else
        re_str = '';
    end
    fname{i} = fullfile(self.procdir,[re_str,self.fn_base,'.',self.dat.(dat_str{i}).tag,'.nii.gz']);
end