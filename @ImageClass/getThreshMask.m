% ImageClass function
% Retrieve thresholded mask
function tmask = getThreshMask(self,vec,dim,slc)
tmask = [];
if self.check && (nargin>1) && ~isempty(vec) && isnumeric(vec)
    vec = round(vec);
    if (nargin == 2) % no slices input, return all
        if self.mask.check
            tmask = self.mask.mat(:,:,:,1);
        else
            tmask = true(self.dims(1:3));
        end
        for i = 1:length(vec)
            tmask = tmask ...
                & (self.mat(:,:,:,vec(i)) >= self.thresh(vec(i),1)) ...
                & (self.mat(:,:,:,vec(i)) <= self.thresh(vec(i),2));
        end
    elseif (nargin==4) && ~isempty(slc) && all(slc>0) ...
            && all(slc<=self.dims(dim)) && (vec>0) && (vec<=self.dims(4))
        slc = round(slc); % make sure integer indices
        tmask = self.mask.getSlice(dim,1,slc);
        for i = 1:length(vec)
            timg = self.getSlice(dim,vec(i),slc);
            tmask = tmask ...
                & (timg >= self.thresh(vec(i),1)) ...
                & (timg <= self.thresh(vec(i),2));
        end
    end
end