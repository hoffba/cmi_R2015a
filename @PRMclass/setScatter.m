% PRMclass function
% Set option to display PRM scatterplot
function setScatter(self,val)
if (nargin == 2) && islogical(val)
    self.prmscatter = val;
end