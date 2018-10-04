function N = excelDate(x)
% Converts Matlab numeric or text date to Excel date number

if ischar(x)
    N = datenum(x,fmt)
    
elseif isnumeric(x)
else
    error('Invalid input type.');
end
