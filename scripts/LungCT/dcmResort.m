% dcmResort.m
% take folder of dcm files and resort into a set of folders with files
% allocated by 'SeriesNumber'
function dcmResort
fnames = dir;
for i = 1:numel(fnames)
    if ~isfolder(fnames(i).name) && isdicom(fnames(i).name)
        info = dicominfo(fnames(i).name);
        serpath = fullfile(pwd,num2str(info.SeriesNumber));
        if ~isfolder(serpath)
            mkdir(serpath)
        end
        movefile(fnames(i).name,serpath,'f');
    else
        continue
    end
end