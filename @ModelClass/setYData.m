% ModelClass function
% Set y-data (Signal Intensity)
function stat = setYData(self,y)
stat = false;
if (nargin==2) && ~isempty(y)
    self.ydata = y(:)';
    %self.updatePlot('data');
    stat = true;
    if length(self.x) ~= length(self.ydata)
        self.x = 0:(length(self.ydata)-1);
    end
end