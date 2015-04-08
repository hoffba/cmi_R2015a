% RegClass function
function setDefVal(self,x,~)
% Set Output directory for all Elastix/Transformix files

if (nargin==1) || ~isempty(x)
    val = str2double(inputdlg('Input default value:','Default Value',1,{num2str(self.defVal)}));
elseif ishandle(x(1)) && strcmp(x.Tag,'edit_defVal')
    val = str2double(x.String);
elseif isnumeric(x)
    val = x;
else
    error('RegClass.setDefVal : Invalid input.');
end

if isempty(val) || isnan(val)
    set(self.h.edit_defVal,'String',num2str(self.defVal));
else
    self.defVal = val;
end
