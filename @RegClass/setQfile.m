% RegClass function
function setQfile(self,x,~)
% Set Output directory for all Elastix/Transformix files

% Validate inputs:
if (nargin==1) || (~isempty(x) && ishandle(x(1)))
    [fname,fpath] = uiputfile('*.txt','Set Queue File:',self.qfile);
    if ischar(fname)
        str = fullfile(fpath,fname);
        if ~strcmp(str(end-3:end),'.txt')
            str = [str,'.txt'];
        end
    else
        str = '';
    end
elseif ischar(x)
    str = x;
else
    warning('Input is not a file path.');
end

% Set RegClass property:
if ~isempty(str)
    self.qfile = str;
    set(self.h.edit_Qfile,'String',str);
end