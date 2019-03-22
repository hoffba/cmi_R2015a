function fname = pathFromFlux(fname)
% Change file name path from local to Flux:

if iscellstr(fname)
    turbo_path = '/nfs/turbo/umms-cgalban/';

    if ischar(fname)
        fname = {fname};
    end
    for i = 1:numel(fname);
        if strcmp(fname{i}(2),':') % PC
            fname{i} = regexprep(fname{i},'\','/');
            fname{i} = strcat(turbo_path,fname{i}(4:end));
        else % Mac
            fname{i} = regexprep(fname{i},'^(.*)/umms-cgalban','/nfs/turbo/umms-cgalban');
        end
    end
else
    error('Invalid input');
end