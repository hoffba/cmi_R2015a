% PRMclass function
% Retrieve PRM region percentages (of total VOI volume)
function [labels,vals] = getStats(self)
if self.check
    np = numel(find(~isnan(self.mat)));
    vals = zeros(1,self.nprm);
    for i = 1:self.nprm
        vals(i) = numel(find(self.mat==i))/np*100;
    end
    labels = self.prmmap(:,2)';
else
    labels = []; vals = [];
end