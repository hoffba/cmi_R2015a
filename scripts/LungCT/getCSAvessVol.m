function [vol] = getCSAvessVol(V, C, csa1, csa2)
%getCSAvessVol.m Creates CSA maps from vessel segmentations
% This is a crude way of getting the vessel segments between a certain CSA
% diameter range. 
% NOTE: Will need to improve this and get the exact computation later
% Inputs:
%   V :    (3D binary image) vessel segmentation
%   C :    CSA Skeleton Image
%   csa1:  CSA value 1 (lower CSA value below which you do not want any vessels)
%   csa2:  CSA value 2 (lower CSA value above which you do not want any vessels)
% Outputs: 
%   vol : Binary Vessels within a certain CSA
    
    %get the CSA skeleton within the range that you want
    C2 = (C<csa2 & C>csa1);
     
    %dilate the range limited CSA skeleton
    % first calculate the diameter of the structing element
    sz = floor(sqrt(csa2/pi)/0.625);
    
    C2dil = imdilate(C2,strel('sphere',sz));
    
    %Get the vessel segments in the required CSA diameter range
    vol = logical(V).* C2dil;
end
%imreconstruct retrieves all connected pieces instead of just the skeleton