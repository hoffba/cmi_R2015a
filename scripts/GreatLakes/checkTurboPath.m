function [fn,flag,LOC_turbo_dir] = checkTurboPath(fnames,turboID)
% Inputs:
%       fnames:     <cellstr> OR <char> files to check
%       turboID:    <char> account ID for Turbo (e.g. 'umms-cgalban') 
% Outputs:
%       fn:         <cellstr> files that pass the test
%       shortpath:  <cell(cellstr)> corresponding path directories split 
%       flag:       <TF> which filenames were not on Turbo

if nargin < 2
    turboID = 'umms-cgalban';
end
GL_turbo_dir = ['/nfs/turbo/',turboID];
LOC_turbo_dir = '';
if ispc
    netdrives = findNetDrives(turboID);
    LOC_turbo_dir = netdrives.Drive;
end

charflag = ischar(fnames);
if charflag
    fnames = {fnames};
end
nf = numel(fnames);
flag = false(nf,1); % flag fnames to remove
fn = cell(nf,1);
for i = 1:nf
    if startsWith(fnames{i},GL_turbo_dir)
        fn{i} = fnames{i};
    else
        str = strsplit(fnames{i},{'\','/'});
        if isempty(LOC_turbo_dir) || strcmp(str{1},LOC_turbo_dir)
            fn{i} = strjoin([GL_turbo_dir,str(2:end)],'/');
        else
            flag(i) = true;
            fprintf('File not on Turbo: %s\n',fnames{i});
        end
    end
end
fn(flag) = [];

if charflag
    fn = fn{1};
end