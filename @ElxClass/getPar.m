% ElxClass function
% Set Elastix parameters in existing schedule
function val = getPar(self,i,fldn)
% Get current parameter value fldn from schedule step i

val = [];
if (nargin==3) && ismember(i,1:length(self.Schedule))
    if ischar(fldn)
        if isfield(self.Schedule{i},fldn)
            val = self.Schedule{i}.(fldn);
        end
    elseif iscellstr(fldn)
        val = cell(1,length(fldn));
        for j = find(isfield(self.Schedule{i},fldn))
            val{j} = self.Schedule{i}.(fldn{j});
        end
    end
end
