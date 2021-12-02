function [fn,shortpath,flag] = checkTurboPath(fnames,dcode)

if ischar(fnames)
    fnames = {fnames};
end
nf = numel(fnames);
flag = false(nf,1); % flag fnames to remove
fn = cell(nf,1);
shortpath = cell(nf,1);
for i = 1:nf
    str = strsplit(fnames{i},'\');
    if ~strcmp(str{1},dcode)
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
