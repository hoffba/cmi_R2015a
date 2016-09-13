% ImageClass function
function setImgOrder(self,i)
if all(ismember(i,1:self.dims(4)))
    
    % Validate image order:
    % remove redundancies
    [~,ia,~] = unique(i);
    i = i(sort(ia));
    % append unmentioned images to back
    i = [i,setdiff(1:self.dims(4),i)];
    
    % Re-order dim4 and associated values:
    self.mat = self.mat(:,:,:,i);
    self.labels = self.labels(i);
    self.valExt = self.valExt(i,:);
    self.scaleM = self.scaleM(i);
    self.scaleB = self.scaleB(i);
    self.thresh = self.thresh(i,:);
    self.fnames = self.fnames(i);
    
else
    error('Requested image number out of range.');
end