% dcmResort.m
% take folder of dcm files and resort into a set of folders with files
% allocated by 'SeriesNumber'
function dcmResort(dcmfolder,destfolder)

if nargin<2
    destfolder = dcmfolder;
end
fn = dir(dcmfolder);
for i = 1:numel(fn)
    fnfull = fullfile(fn(i).folder,fn(i).name);
    if ~fn(i).isdir && isdicom(fnfull)
        info = dicominfo(fnfull);
        serpath = fullfile(destfolder,sprintf('ser%d',info.SeriesNumber));
        if ~isfolder(serpath)
            mkdir(serpath)
        end
        try
            movefile(fnfull,serpath,'f');
        catch err
            disp('break')
        end
    end
end