function [seg,np,T] = quantAirTrapping(ins,imask,exp,emask)
% Method for quantifying air trapping (Goris et al. 2003)
% Segments lungs based on median, 90%ile of histogram
% Inputs:
%       ins = inspiratory image
%       imask = inspiratory binary mask
%       exp = expiratory image
%       emask = expiratory binary mask
% Output:
%       seg = vector of indexes for segmentation
%       np = segmentation number of voxels

% Default HU limit:
HUmax = -250;

% Determine thresholds
ipctl = prctile(ins(imask & (ins<HUmax)),[90,50]);
epctl = prctile(exp(emask & (exp<HUmax)),90);
D = epctl-ipctl(1);
T = ipctl(1) - ((0:2) - (1-D/343)) * diff(ipctl)/3;

% Segmentation of the expiratory image:
ns = size(exp,3);
np = zeros(1,3);
seg = 0*exp;
hw = waitbar(0,'Segmenting expiratory image ...');
for i = 1:3
    tseg = exp < T(i);
    for j = 1:ns
        tseg(:,:,j) = medfilt2(tseg(:,:,j),[15,15]);
        waitbar(j/ns,hw,['Segmenting expiratory image ... ',num2str(i)]);
    end
    np(i) = nnz(tseg);
    seg = seg + tseg;
end
delete(hw);
