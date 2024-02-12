% Raw 3D MF calculation
function p = calcMF3D(BW,voxsz,ord,mask)
if nargin==4
    normval = imVolumeEstimate(mask,voxsz);
else
    normval = 1;
end
p = nan(1,length(ord));
if ismember(0,ord)
    p(ord==0) = imVolumeEstimate(BW,voxsz)/normval;
end
if ismember(1,ord)
    p(ord==1) = imSurfaceEstimate(BW,voxsz)/normval;
end
if ismember(2,ord)
    p(ord==2) = imMeanBreadth(BW,voxsz)/normval;
end
if ismember(3,ord)
    p(ord==3) = imEuler3dEstimate(BW)/nnz(mask);
end
