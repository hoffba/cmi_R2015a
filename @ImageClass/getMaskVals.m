% ImageClass function
function vals = getMaskVals(self,vec)
vals = [];
if (nargin==1)
    vec = 1:self.dims(4);
end
if self.check && self.mask.check && all(vec>0) && all(vec<=self.dims(4))
    np = sum(self.mask.mat(:));
    vals = zeros(np,length(vec));
    for i = 1:length(vec)
        timg = self.mat(:,:,:,vec(i));
        vals(:,i) = timg(self.mask.mat);
    end
end