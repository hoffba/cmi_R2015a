% HistogramClass function
% Set active status
function setActive(self,val)
self.active = logical(val);
self.dispUDvisible;