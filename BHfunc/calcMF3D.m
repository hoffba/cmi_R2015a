% Raw 3D MF calculation
function p = calcMF3D(BW,voxsz,ord,mask)
if nargin==4
    normvol = imVolumeEstimate(mask,voxsz);
    normN = nnz(mask);
else
    normvol = 1;
    normN = 1;
end
p = nan(1,length(ord));
if ismember(0,ord)
    p(ord==0) = imVolumeEstimate(BW,voxsz)/normvol;
end
if ismember(1,ord)
    p(ord==1) = imSurfaceEstimate(BW,voxsz)/normvol;
end
if ismember(2,ord)
    p(ord==2) = imMeanBreadth(BW,voxsz)/normvol;
end
if ismember(3,ord)
    p(ord==3) = imEuler3dEstimate(BW)/normN;
end
