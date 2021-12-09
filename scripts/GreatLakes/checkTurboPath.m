function [fn,shortpath,flag,drivebase] = checkTurboPath(fnames,drivebase)
% Inputs:
%       fnames:     <cellstr> OR <char> files to check
%       drivebase:  <char> local path to Turbo 
% Outputs:
%       fn:         <cellstr> files that pass the test
%       shortpath:  <cell(cellstr)> corresponding path directories split 
%       flag:       <TF> which filenames were not on Turbo

if nargin < 2
    netdrives = findNetDrives('umms-cgalban');
    drivebase = netdrives.Drive;
end

if ischar(fnames)
    fnames = {fnames};
end
nf = numel(fnames);
flag = false(nf,1); % flag fnames to remove
fn = cell(nf,1);
shortpath = cell(nf,1);
for i = 1:nf
    str = strsplit(fnames{i},filesep);
    if ~strcmp(str{1},drivebase)
        flag(i) = true;
        fprintf('File not on Turbo: %s\n',fnames{i});
    else
        if isfolder(fnames{i})
            if nf == 1
                shortpath{i} = str;
            else
                flag(i) = true;
                fprintf('Value must be a file path, not directory: %s\n',fnames{i});
            end
        else
            shortpath{i} = str(1:end-1);
            fn(i) = str(end);
        end
    end
end
fn(flag) = [];
shortpath(flag) = [];
