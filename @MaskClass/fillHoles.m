% MaskClass function
% Fill holes in current VOI
function fillHoles(self)
if self.check
    tdims = self.dims;
    for i = 1:tdims(3)
        self.mat(:,:,i,1) = imfill(self.mat(:,:,i,1),'holes');
    end
end