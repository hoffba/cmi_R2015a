% ImageClass function
% Delete Image
function imgDelete(self,inds)
% inds can be a vector of indices within range of the image matrix
% vec is optional for keeping track of a certain image
if (nargin > 1)
    if isnumeric(inds) && ~isempty(inds)
        inds = unique(uint8(inds));
        nvec = self.dims(4);
        if ~any(inds<1 | inds>nvec)
            self.mat(:,:,:,inds) = [];
            self.dims(4) = nvec - length(inds);
            self.labels(inds) = [];
            self.thresh(inds,:) = [];
            self.scaleM(inds) = [];
            self.scaleB(inds) = [];
            self.valExt(inds,:) = [];
        end
    end
end