% ImageClass function
% Retrieve thresholded mask
function tmask = getThreshMask(self,vec,dim,slc)
tmask = [];
if self.check && (nargin>1) && ~isempty(vec) && isnumeric(vec)
    vec = round(vec(1));
    if (nargin == 2) % no slices input, return all
        if self.mask.check
            tmask = self.mask.mat(:,:,:,1);
        else
            tmask = ones(self.dims(1:3));
        end
        timg = self.mat(:,:,:,vec);
    elseif (nargin==4) && ~isempty(slc) && all(slc>0) ...
            && all(slc<=self.dims(dim)) && (vec>0) && (vec<=self.dims(4))
        slc = round(slc); % make sure integer indices
        tmask = self.mask.getSlice(dim,1,slc);
        timg = self.getSlice(dim,vec,slc);
    end
    if ~isempty(tmask)
        tmask = tmask ...
              & (timg >= self.thresh(vec,1)) ...
              & (timg <= self.thresh(vec,2));
    end
end