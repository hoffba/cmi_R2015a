% Raw 3D MF calculation
function p = calcMF3D(BW,voxsz,ord,mask)
if nargin==4
    maskvol = imVolumeEstimate(mask,voxsz);
    masknum = numel(mask);
else
    maskvol = 1;
    masknum = 1;
end
p = nan(1,length(ord));
if ismember(0,ord)
    p(ord==0) = imVolumeEstimate(BW,voxsz)/maskvol;
end
if ismember(1,ord)
    p(ord==1) = imSurfaceEstimate(BW,voxsz)/maskvol;
end
if ismember(2,ord)
    p(ord==2) = imMeanBreadth(BW,voxsz)/maskvol;
end
if ismember(3,ord)
    p(ord==3) = imEuler3dEstimate(BW)/masknum;
end
