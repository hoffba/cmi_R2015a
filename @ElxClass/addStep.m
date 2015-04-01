% ElxClass function
% Set Initial Transform Options
function stat = addStep(self,str,varargin)

stat = false;
if nargin<2
    str = 'affine';
end
i = length(self.Schedule)+1;
p = initElastixParameters(str);
if ~isempty(p)
    self.Schedule{i} = p;
    self.setPar(i,varargin{:});
    stat = true;
end