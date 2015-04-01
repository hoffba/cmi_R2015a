% ImageClass function
function [vMean,vStD,vMed,vVol,nVox] = getStats(self,vec)
% Calculates the following VOI statistics:
%       Mean Value
%       Standard Deviation
%       Median Value
%       Volume
%       # Voxels

vMean = [];
vStD = [];
vMed = [];
vVol = [];
nVox = [];
% only perform stats on loaded image using current mask
if self.check && self.mask.check
    if nargin<2
        vec = 1:self.dims(4);
    elseif isnumeric(vec)
        vec = round(vec);
        vec = vec((vec>0) & (vec<=self.dims(4)));
    else
        vec = [];
    end
    nv = length(vec);
    vMean = zeros(1,nv); vStD = vMean; vMed = vMean; vVol = vMean; nVox = vMean;
    for i = 1:length(vec)
        tvals = self.mat(:,:,:,vec(i));
        tvals = tvals(self.mask.mat);
        vMean(i) = mean(tvals);
        vStD(i) = std(tvals);
        vMed(i) = median(tvals);
        nVox(i) = numel(tvals);
        vVol(i) = nVox(i) * prod(self.voxsz);
    end
end