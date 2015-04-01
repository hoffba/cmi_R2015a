function [files,paths] = cmidir(filt,sub,fcheck)
% directory search for files
% inputs:    filt = sting filter for file search
%            sub = option to search sub-directories
%            fcheck = option for GUI to select files
% output:    cell array of file names

fpath = uigetdir('Select folder to search:');
if fpath
    if sub == 1
        t = subsearch(fpath,filt);
        if isempty(t)
            error('No files were found matching input criteria.');
        else
            files = cell(size(t,1),1); paths = files;
            for i = 1:size(t,1)
                [folder, name, ext] = fileparts(t(i).name);
                files{i} = [name ext];
                paths{i} = [folder filesep];
            end            
        end
    else
        cd(fpath)
        t = dir(filt);
        if isempty(t)
            error('No files were found matching input criteria.');
        else
            files = cell(size(t,1),1); paths = files;
            for i = 1:size(t,1)
                files{i} = t(i).name;
                paths{i} = [fpath filesep];
            end
        end
    end
else
    files = []; paths = [];
end
if fcheck % provides GUI for selecting files
    [sel,ok] = listdlg('PromptString','Select files to process:',...
        'ListSize',[250 300],'ListString',files);
    if ok
        files = files(sel);
        paths = paths(sel);
    else
        files = {}; paths = {};
    end
end

function Files = subsearch(folder,filter)
% searches all subfolders
pathstr = genpath(folder);
seplocs = strfind(pathstr, pathsep);
loc1 = [1 seplocs(1:end-1)+1];
loc2 = seplocs(1:end)-1;
pathfolders = arrayfun(@(a,b) pathstr(a:b), loc1, loc2, 'UniformOutput', false);
Files = [];
for ifolder = 1:length(pathfolders)
    NewFiles = dir(fullfile(pathfolders{ifolder}, filter));
    if ~isempty(NewFiles)
        fullnames = cellfun(@(a) fullfile(pathfolders{ifolder}, a), {NewFiles.name},...
            'UniformOutput', false); 
        [NewFiles.name] = deal(fullnames{:});
        Files = [Files; NewFiles];
    end
end