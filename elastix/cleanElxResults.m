function cleanElxResults(fdir)
% Deletes result.* files from directory tree

if nargin==0
    fdir = uigetdir(pwd,'Select parent directory for result.* search:');
end
if ischar(fdir) && isdir(fdir)
    [D,F] = dirtree(fdir,'result.*');
    if ~isempty(D)
        fnames = cellfun(@(x,y)fullfile(x,y)',D,F,'UniformOutput',false);
        fnames = [fnames{:}]';
        hw = waitbar(0,'Deleting files ...');
        nf = length(fnames);
        for i = 1:nf
            try
                delete(fnames{i})
                fnames{i} = '';
            catch
            end
            waitbar(i/nf,hw);
        end
        delete(hw);
        fnames(cellfun(@isempty,fnames)) = [];
        if ~isempty(fnames)
            fprintf('Could not delete:\n');
            fprintf('  %s\n',fnames{:});
        end
    end
end