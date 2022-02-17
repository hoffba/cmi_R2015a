% Raw 2D MF calculation
function p = calcMF2D(BW,voxsz,ord)
p = nan(1,length(ord));
if ismember(0,ord)
    p(ord==0) = imAreaEstimate(BW,voxsz);
end
if ismember(1,ord)
    p(ord==1) = imPerimeterEstimate(BW,voxsz);
end
if ismember(2,ord)
    p(ord==2) = imEuler2dEstimate(BW);
end