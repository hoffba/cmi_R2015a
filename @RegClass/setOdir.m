% RegClass function
function setOdir(self,x,~)
% Set Output directory for all Elastix/Transformix files

if nargin==1
    x = self.h.button_odir;
end

if ischar(x)
    str = x;
elseif ishandle(x(1))
    switch x(1).Tag
        case 'button_odir'
            odir = uigetdir(self.odir,'Select output directory');
            if ischar(odir)
                str = odir;
            else
                str = '';
            end
        case 'edit_odir'
            str = x.String;
        otherwise
            str = '';
    end
else
    warning('RegClass/setOdir: Invalid input.');
end

if ~isempty(str)
    % Check for ouput directory:
    if self.guicheck
%         if ~exist(str,'dir')
%             answer = questdlg(['Directory does not exist. Create directory: ',str,'?']);
%             if strcmp(answer,'Yes')
%                 mkdir(str);
%             else
%                 str = '';
%             end
%         end
        set(self.h.edit_Odir,'String',str);
    end
    self.odir = str;
end