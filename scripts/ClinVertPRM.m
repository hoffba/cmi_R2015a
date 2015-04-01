% CMI script
function ClinVertPRM(cmiObj)
% Performs PRM analysis on selected data: from first scan and sequentailly

%   Input: cmiObj = CMIclass object containing current settings

% Select images to process
fnames = {};
tdir = pwd;
while tdir~=0
    tdir = uigetdir(tdir,'Select Data Folder:');
    if ischar(tdir)
        fnames{end+1} = tdir;
        [] = fileparts(tdir);
    end
end