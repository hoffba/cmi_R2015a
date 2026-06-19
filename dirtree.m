function [D,F] = dirtree(tpath,filtstr)

% Loops over filters
if ischar(filtstr)
    filtstr = {filtstr};
end
fn = dir(fullfile(tpath,'**'));

% Remove directories
fn([fn.isdir]) = [];

% Remove files not matching filtstr
ind = endsWith({fn.name},filtstr);
if ismember(filtstr,'')
    ind = ind | ~contains({fn.name},'.');
end
fn(~ind) = [];

% Return directory (D) and filename (F) cell arrays
[D,~,ic] = unique({fn.folder}');
ndir = numel(D);
F = cell(ndir,1);
for i = 1:numel(D)
    f = {fn(ic==i).name}';
    F{i} = unique(f);
end

% Remove search folder if there are subfolders
if numel(D)>1
    ind = strcmp(D,tpath);
    D(ind) = [];
    F(ind) = [];
end