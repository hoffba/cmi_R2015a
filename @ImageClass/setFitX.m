% ImageClass function
function stat = setFitX(self,x)
% Set XData for curve-fitting
% OUTPUT: stat - success(true), fail(false), cancel(empty)
stat = false;
if self.check
    if nargin==1 % GUI input
        x = [];
        if isempty(self.model.xdata)
            def = cellfun(@str2double,self.labels);
        else
            def = num2str(self.model.xdata);
        end
        if any(isnan(def))
            def = {''};
        else
            def = {num2str(def)};
        end
        answer = inputdlg('Input X-values (incr/vector):','XData',1,def);
        if isempty(answer)
            stat = [];
        else
            try
                x = eval(answer{1});
            catch err
                rethrow(err);
            end
        end
    end
    nx = length(x);
    if any(nx==[1,self.dims(4)]) && isnumeric(x) && ~any(isnan(x))
        if nx==1
            x = x*(0:(self.dims(4)-1));
        end
        stat = self.model.setXData(x);
    end
end