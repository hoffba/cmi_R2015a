% ElxClass function
% Set Elastix parameters in existing schedule
function stat = setPar(self,i,varargin)
% Set current parameter value
% OR add new parameter value
% OR remove existing parameter (value of Name/Value pair is empty)

stat = false;

if (nargin>2) && ismember(i,1:length(self.Schedule))
    if isstruct(varargin{1})
        tstruct = varargin{1};
        str = fieldnames(tstruct);
        for j = 1:length(str)
            if isfield(self.Schedule{i},str{j})
                self.Schedule{i}.(str{j}) = tstruc.(str{j});
            end
        end
    elseif mod(length(varargin),2)==0 % Name/Value pairs
        for j = 1:2:length(varargin)
            if ischar(varargin{j})
                if isempty(varargin{j+1}) && isfield(self.Schedule{i},varargin{j})
                    self.Schedule{i} = rmfield(self.Schedule{i},varargin{j});
                else
                    self.Schedule{i}.(varargin{j}) = varargin{j+1};
                end
            end
        end
    else
        error('Invalid inputs.')
    end
end
