% RegClass function
function setQnum(self,hObject,~)

% Parse input:
if isa(hObject,'matlab.ui.control.UIControl')
    ival = str2double(hObject.String);
elseif isnumeric(hObject)
    ival = hObject;
else
    ival = nan;
end

% Set Qnum value:
nmax = feature('NumCores');
ival = round(ival);
if ~(isnan(ival) || isempty(ival)) && (ival>0)
    if ival>nmax
        warning(['Number of jobs set to maximum: ',num2str(nmax)]);
        ival = nmax;
    end
    self.qnum = ival;
else
    self.h.edit_Qnum.String = num2str(self.qnum);
    warning('Invalid input: number of jobs.');
end
self.h.edit_Qnum.String = num2str(self.qnum);
