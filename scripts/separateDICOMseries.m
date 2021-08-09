
function separateDICOMseries(bpath)

% Separates DICOM files into series folders

if nargin==0
    bpath = pwd;
end

dnames = dir(bpath);
dnames = {dnames(cellfun(@isfolder,{dnames.name})).name}';
dnames(1:2) = [];
answer = listdlg('ListString',dnames);
dnames = dnames(answer);
for idir = 1:length(dnames)

    tpath = fullfile(bpath,dnames{idir});
    fnames = dir(fullfile(tpath,'*.dcm'));
    for i = 1:length(fnames)
        info = dicominfo(fullfile(tpath,fnames(i).name));
        serUID = info.SeriesInstanceUID;
        if ~isfolder(fullfile(tpath,serUID))
            mkdir(fullfile(tpath,serUID));
        end
        movefile(fullfile(tpath,fnames(i).name),fullfile(tpath,serUID));
    end
    
end