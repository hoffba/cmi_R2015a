function newp = cmi_setPath(path_cmi)
% Compiles list of directories to add to the Matlab path for the CMI
% program
% 1. Run this script from the cmi_R2015a directory

if nargin==0
    path_cmi = pwd;
end
if ~strcmp(path_cmi(end-9:end),'cmi_R2015a')
    error('Must run this script from the cmi_R2015a directory.')
end
newp = genpath(pwd);
newp = strsplit(newp,pathsep);

% Remove .git folders:
newp(contains(newp,'.git')) = [];

% Remove OLD folders:
newp(contains(newp,'OLD')) = [];

% Remove directories already existing in the search path
p = strsplit(path,pathsep);
newp(ismember(newp,p)) = [];

% Add to Matlab path
addpath(newp{:});

fprintf('\n\nAdded folders to Matlab search path:\n')
fprintf('   %s\n',newp{:});
fprintf('\n\n');