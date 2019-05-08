% ImageClass function
% Set labels for 4D images
function setLabel(self,vec,str)
if ischar(str)
    if numel(vec)==1
        str = {str};
    else
        error('Invalid inputs.')
    end
end
vec = round(vec);
if length(vec)~=length(str)
    error('Number of vectors must equal number of label strings.');
elseif self.check && isnumeric(vec) && iscellstr(str)
    self.labels(round(vec)) = str;
end