% ImageIO class function
% [M,info] = read()
% [M,info] = read(fname)
% [M,info] = read(fname,'Name',Value)
% [M,info] = read('Name',Value)
% 
% Name/Value pairs:
%   'd' : matrix dimensions to match

function [M,cmiInfo,natInfo] = read(self,varargin)

M = [];
cmiInfo = [];
natInfo = [];

n = length(varargin);
if (n==0) || ~exist(varargin{1},'file')
    fname = {''};
else
    fname = varargin{1};
    varargin(1) = [];
end

% Determine inputs:
p = inputParser;
addRequired(p,'fname',@(x)ischar(x)||iscellstr(x));
addParameter(p,'d',[],@isvector);
addParameter(p,'filt',strcat('*',[self.types(:).ext]'),@iscellstr);
parse(p,fname,varargin{:});
p = p.Results;
if ischar(p.fname)
    p.fname = {p.fname};
end

if self.gui
    % Use input name as default for selecting image:
    if iscellstr(p.fname)
        p.fname = p.fname{1};
    end
    [fname,fpath] = uigetfile(self.getFilt,'Load Image:',p.fname,...
        'MultiSelect','on');
    if ischar(fpath)
        p.fname = fullfile(fpath,fname);
    else
        p.fname = {};
    end
end

if ischar(p.fname)
    p.fname = cellstr(p.fname);
end
nf = length(p.fname);
for i = 1:nf
    % Find load function:
    [func,itype] = self.getFunc(p.fname{i});
    if isempty(func)
        disp('Load function not yet programmed.');
    else
        % Load image:
        [M,info] = feval(func,p.fname{i},p.d);
        if isempty(M)
            return;
        end
        if i==1
            % Only save info for first one loaded
            natInfo = info;
            cmiInfo = self.parseInfo(info,itype);
            if isempty(p.d)
                p.d = cmiInfo.Dim(1:3);
            end
        end
        if i==nf
            % Remember format of last image loaded:
            self.wdir = fileparts(p.fname{i});
            self.lasttype = itype;
        end
    end
end