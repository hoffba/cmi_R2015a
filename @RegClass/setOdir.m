% RegClass function
function setOdir(self,x,~)
% Set Output directory for all Elastix/Transformix files

if (nargin==1) || (~isempty(x) && ishandle(x(1)))
    odir = uigetdir(self.odir,'Select output directory');
    if ischar(odir)
        str = odir;
    else
        str = '';
    end
elseif ischar(x)
    str = x;
else
    warning('Input is not a directory path.');
end

if ~isempty(str)
    % Check for ouput directory:
    if ~exist(str,'dir')
        answer = questdlg(['Directory does not exist. Create directory: ',str,'?']);
        if strcmp(answer,'Yes')
            mkdir(str);
        else
            str = '';
        end
    end
    if ~isempty(str)
        self.odir = str;
        set(self.h.edit_Odir,'String',str);
    end
end