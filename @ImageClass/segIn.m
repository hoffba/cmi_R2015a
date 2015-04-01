% ImageClass function
% Grow VOI from one ROI to another using threshold value
function segIn(self,vec,step)
% Requires that an ROI be drawn on each end of what you want segmented
% Will fill in the middle through dilation / threshold
% Inputs: thresh = threshold value for segmenting desired region [min,max]
%         step   = dilation radius [x-y,z]
if ~isempty(self.mat) && ~isempty(self.mask.mat) && (nargin == 3) ...
        && isnumeric(vec) && isnumeric(step) && (length(step)==2)
    tslc = find(squeeze(any(any(self.mask.mat)))); % To find ends
    if length(tslc) > 1
        tslc = [min(tslc)-1 max(tslc)+1];
        thresh = self.thresh(vec,:);
        % thmask defines all voxels within the slice range and within the threshold range
        thmask = (self.mat(:,:,:,vec) >= thresh(1)) & (self.mat(:,:,:,vec) <= thresh(2));
        thmask(:,:,[1:tslc(1),tslc(2):end]) = false;
        self.mask.merge('intersect',thmask); % threshold original mask
        done = false;
        count = 1;
        while ~done
            disp(['step ' num2str(count)])
            omask = self.mask.mat;
            self.mask.morph('dilate',step); % grow the mask
            self.mask.merge('intersect',thmask); % threshold the mask
            if all(omask(:) == self.mask.mat(:)) % check that mask has changed
                done = true;
            end
            count = count + 1;
        end
    end
end