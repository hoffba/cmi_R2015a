% ElxClass function
% Set Elastix parameters in existing schedule
function stat = setPar(self,i,varargin)
% Set current parameter value
% OR add new parameter value
% OR remove existing parameter (value of Name/Value pair is empty)

stat = false;

if (nargin>2) && all(ismember(i,1:length(self.Schedule)))
    if isstruct(varargin{1})
        tstruct = varargin{1};
        str = fieldnames(tstruct);
        for j = 1:length(str)
            if isfield(self.Schedule{i},str{j})
                for ii = 1:length(i)
                    self.Schedule{i(ii)}.(str{j}) = tstruc.(str{j});
                end
            end
        end
        stat = true;
    elseif mod(length(varargin),2)==0 % Name/Value pairs
        for j = 1:2:length(varargin)
            fstr = varargin{j};
            fval = varargin{j+1};
            if ischar(fstr)
                for ii = 1:length(i)
                    ind = i(ii);
                    if strncmp(fstr,'FinalGridSpacingIn',18)
                        tstr = {'Voxels','PhysicalUnits'};
                        tstr = tstr{strcmp(fstr(19:end),tstr)};
                        if isempty(tstr)
                            warning('Unknown tag: %s',fstr);
                        elseif isfield(self.Schedule{ind},tstr)
                            self.Schedule{ind} = rmfield(self.Schedule{ind},tstr);
                        end
                    end
                    if isempty(fval) && isfield(self.Schedule{ind},fstr)
                        self.Schedule{ind} = rmfield(self.Schedule{ind},fstr);
                    else
                        self.Schedule{ind}.(fstr) = fval;
                    end
                end
            end
        end
        stat = true;
    else
        error('Invalid inputs.')
    end
end
