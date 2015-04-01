% ImageClass function
% Add threshold to mask
function thresh2mask(self,vec)
tmask = self.getThreshMask(vec); % only one input for 3D mask
if ~isempty(tmask)
    if self.mask.check
        str = 'intersect';
    else
        str = 'add';
    end
    self.mask.merge(str,tmask);
    self.threshold(vec,10^6*[-1 1]');
end