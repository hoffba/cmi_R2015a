% ElxClass function
% Remove schedule step
function stat = rmStep(self,i)

stat = true;
if nargin==1
    % Remove all
    self.Schedule = {};
elseif (nargin==2) && ismember(i,1:length(self.Schedule))
    self.Schedule(i) = [];
else
    stat = false;
end