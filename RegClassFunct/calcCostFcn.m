% Calculates cost function for image coregistration
% e = calcCostFcn(x,y,metric)
% Inputes:  x = vector of reference values
%           y = vector of homologous values
%           metric = cost function option
%                   mi : Mutual Information
%                   sd : Squared Difference
function e = calcCostFcn(x,y,metric)
e = [];
if (nargin>1) && (length(x)==length(y))
    if ~exist('metric','var') || ~ischar(metric)
        metric = sd;
    end
    nx = length(x);
    switch metric
        case mi
            
        case sd
            e = sum((x-y).^2)/nx;
    end
end
