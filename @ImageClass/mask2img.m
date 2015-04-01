% ImageClass function
function stat = mask2img(self,str)
% Appends current mask to end of 4D image

stat = false;
if self.mask.check
    
    % optional manual input str = image label
    if ~ischar(str)
        str = 'VOI';
    end
    
    self.mat(:,:,:,end+1) = double(self.mask.mat);
    self.labels{end+1} = str;
    self.scaleM(end+1) = 1;
    self.scaleB(end+1) = 0;
    self.thresh(end+1,:) = [-inf,inf];
    self.dims(4) = size(self.mat,4);
    self.valExt(end+1,:) = [0,1];
    
    self.mask.clear(0); % clear the mask
    
    stat = true;
end


