function [D,F] = dirtree(tpath,filtstr)

% Loops over filters
if ischar(filtstr)
    filtstr = {filtstr};
end
fn = [];
for i = 1:numel(filtstr)
    fn = [fn;dir(fullfile(tpath,['**',filesep,filtstr{i}]))]; %#ok<AGROW>
end

% Remove directories
fn([fn.isdir]) = [];

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