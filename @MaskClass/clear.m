% MaskClass function
% Clear mask (slice/all)
function clear(self,slc,sview)

if (nargin==1)
    self.mat = [];
    self.check = false;
elseif (nargin==3) && (slc>0) && any(sview==1:3) && (slc<=self.dims(sview))
    tdims = self.dims(1:3);
    tdims(sview) = [];
    nslc = zeros(tdims);
    self.setSlice(nslc,sview,1,slc)
    if ~any(self.mat(:))
        self.mat = [];
        self.check = false;
    end
end
