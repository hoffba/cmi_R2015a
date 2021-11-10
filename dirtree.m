function [D,F] = dirtree(tpath,filt,h)
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
if (nargin<3) || ~ishandle(h)
    h = false;
%     h = waitbar(0,'Finding directories');
end

D = {};
F = {};
if ~isempty(filt)
    [~,dname] = fileparts(tpath);
    disp(dname);
%     waitbar(0,h,dname);
    tD = dir(tpath);
    % Match file filters in this path:
        fchk = false;
        fnames = {};
        for i = 1:length(filt)
            tfnames = dir(fullfile(tpath,filt{i}));
            % Ignore directories
            tfnames([tfnames(:).isdir]) = [];
            tfnames = {tfnames(:).name};
            % Ignore DICOMDIR
            tfnames(cellfun(@(x)strcmp(x,'DICOMDIR'),tfnames)) = [];
            if isempty(filt{i})
                % Ignore files with extensions
                tfnames(cellfun(@(x)ismember('.',x),tfnames)) = [];
            end
            if ~isempty(tfnames)
                fnames = [fnames;tfnames'];
                fchk = true;
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
            [DD,FF] = dirtree(fullfile(tpath,tD{i}),filt,h);
            D = [ D ; DD];
            F = [ F ; FF];
        end
    end
end

if (nargin<3) && ishandle(h)
%     delete(h);
end