% ModelClass function
% Set x-data (DCE: time, Diff: b-value)
function stat = setXData(self,x)
stat = false;
if nargin==1 % use GUI to input x-values
    x = [];
    answer = inputdlg('Input X-Data:','X-Data',1);
    if ~isempty(answer)
        try
            x = eval(answer);
        catch err
            rethrow(err);
        end
    end
end
nx = length(x);
if (nx>1) && isnumeric(x) && ~any(isnan(x))
    self.x = x;
    stat = true;
    % Make sure ydata is same length
    if length(self.ydata) ~= nx
        self.ydata = zeros(1,nx);
    end
end