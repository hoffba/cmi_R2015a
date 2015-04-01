% CMIclass function
% Set apply all check
function setApplyAllCheck(self,x,~)
val = [];
if (nargin>1)
    if (nargin==3) && ishandle(x) && strcmp(get(x,'Tag'),'applytoall')
        val = logical(get(x,'Value'));
    elseif ~isempty(x) && isnumeric(x)
        val = logical(x);
    end
end
if ~isempty(val)
    self.applyallcheck = val;
    if val % Set all 4D color limits to current image's values
        self.clim(:,1) = self.clim(self.vec,1);
        self.clim(:,2) = self.clim(self.vec,2);
    end
end
