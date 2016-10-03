% ImageClass method
function stat = reorder(self,v)
% Reorder 4D

stat = false;
v = unique(v,'stable'); % don't allow repetitions but keep in order
if all(ismember(v,1:self.dims(4)))
    % Find unmentioned images and add to end:
    v = [v,setdiff(1:self.dims(4),v)];
    self.mat = self.mat(:,:,:,v);
    self.scaleM = self.scaleM(v);
    self.scaleB = self.scaleB(v);
    self.thresh = self.thresh(v,:);
    self.labels = self.labels(v);
    stat = true;
else
    warning('Input indices are outside image limits.');
end