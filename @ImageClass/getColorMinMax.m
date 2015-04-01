% ImageClass function
% Retrieve intensity min/max for current image
function [vmin,vmax] = getColorMinMax(self,vec)
if nargin == 1
    vec = 1:self.dims(4); % return all
elseif all(vec >= 1) && all(vec <= self.dims(4))
    vec = int8(vec);
else
    vec = [];
end
if ~isempty(vec)
    vmin = self.valExt(vec,1);
    vmax = self.valExt(vec,2);
end