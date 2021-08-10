
function separateDICOMseries(bpath)

% Separates DICOM files into series folders

if nargin==0
    bpath = pwd;
end

% Find subdirectories to the base path
n = length(bpath)+1;
dnames = strsplit(genpath(bpath),';')';
dnames = cellfun(@(x)x(n+1:end),dnames,'UniformOutput',false);
dnames(cellfun(@isempty,dnames)) = [];
answer = listdlg('ListString',dnames);
dnames = dnames(answer);

% Perform separation on files in each folder:
for idir = 1:length(dnames)

    tpath = fullfile(bpath,dnames{idir});
    
    % Find DICOM files:
    fnames = dir(fullfile(tpath,'*.dcm'));
    if isempty(fnames)
        % Check for files with no extension
        fnames = dir(tpath);
        fnames(contains({fnames.name},'.')) = [];
        fnames([fnames.isdir]) = [];
        if isempty(fnames)
            fprintf('No DICOM files found in directory: %s\n',dnames{idir});
        end
    end
    
    % Move DICOM files to new subfolders
    for i = 1:length(fnames)
        shortname = fullfile(dnames{idir},fnames(i).name);
        if ~strcmp(fnames(i).name,'DICOMDIR')
            fullname = fullfile(tpath,fnames(i).name);
            if isdicom(fullname)
                info = dicominfo(fullname);
                svdir = info.SeriesNumber;
                if ~isfolder(fullfile(tpath,svdir))
                    mkdir(fullfile(tpath,svdir));
                end
                fprintf('Moving file: %s\n',shortname);
                movefile(fullname,fullfile(tpath,svdir));
            else
                fprintf('Not a DICOM: %s\n',shortname);
            end
        end
    end
    
end