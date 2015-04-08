% RegClass function
function setT0DefVal(self,x,~)
% Set Output directory for all Elastix/Transformix files

if (nargin==1) || isempty(x)
    val = str2double(inputdlg('Input default value:','Default Value',1,{num2str(self.T0defVal)}));
elseif ishandle(x(1)) && strcmp(x.Tag,'edit_T0defVal')
    val = str2double(x.String);
elseif isnumeric(x)
    val = x;
else
    error('RegClass.setDefVal : Invalid input.');
end

if isempty(val) || isnan(val)
    set(self.h.edit_T0defVal,'String',num2str(self.T0defVal));
else
    self.T0defVal = val;
end
