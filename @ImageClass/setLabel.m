% ImageClass function
% Set labels for 4D images
function setLabel(self,vec,str)
if self.check && isnumeric(vec) && ischar(str)
    self.labels{round(vec(1))} = str;
end