% ElxClass function
% Remove schedule step
function stat = rmStep(self,i)

stat = false;
if (nargin==2) && ismember(i,1:length(self.Schedule))
    self.Schedule(i) = [];
    stat = true;
end