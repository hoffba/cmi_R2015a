% dcmResort.m
% take folder of dcm files and resort into a set of folders with files
% allocated by 'SeriesNumber'
function dcmResort(dcmfolder)
fn = dir(dcmfolder);
for i = 1:numel(fn)
    fnfull = fullfile(fn(i).folder,fn(i).name);
    if ~fn(i).isdir && isdicom(fnfull)
        info = dicominfo(fnfull);
        serpath = fullfile(dcmfolder,num2str(info.SeriesNumber));
        if ~isfolder(serpath)
            mkdir(serpath)
        end
        movefile(fnfull,serpath,'f');
    end
end