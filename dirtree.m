function [D,F] = dirtree(tpath,filt)
% Returns full directory tree from base directory tpath of folders 
%   containing files within the filtered set
% Inputs:
%       tpath : string, base directory for search
%       filt : string or cell array of strings, 
%                   filename filter (using wildcards)
%                   Note: Empty string ('') looks for files without any extension,
%                         but empty cell array looks for all directories
% Output: 
%       D : cell array of strings containing directories
%       F : cell array of cell array of strings containing files matching filters

if nargin<2
    filt = {};
elseif ischar(filt)
    filt = {filt};
elseif ~iscellstr(filt)
    filt = {};
end
if nargin==0
    tpath = pwd;
end

D = {};
F = {};
if ~isempty(filt)
%     disp(tpath)
    tD = dir(tpath);
    % Match file filters in this path:
    if isempty(filt)
        fchk = true;
    else
        fchk = false;
        fnames = {};
        for i = 1:length(filt)
            tfnames = dir(fullfile(tpath,filt{i}));
            if isempty(filt{i})
                % Ignore '.' and '..'
                tfnames(1:2) = [];
                % Ignore directories
                tfnames([tfnames(:).isdir]) = [];
                % Ignore files with extensions
                tfnames(cellfun(@(x)ismember('.',x),{tfnames(:).name})) = [];
            end
            fnames = [fnames;{tfnames(:).name}'];
            fchk = fchk || ~isempty(tfnames);
        end
    end
    if fchk
        D = {tpath};
        F = {fnames};
    end
    % Find directories in this path:
    tD = {tD([tD(:).isdir]).name}';
    tD(1:2) = []; % remove "." and ".."
    if ~isempty(tD)
        for i = 1:length(tD)
            [DD,FF] = dirtree(fullfile(tpath,tD{i}),filt);
            D = [ D ; DD];
            F = [ F ; FF];
        end
    end
end