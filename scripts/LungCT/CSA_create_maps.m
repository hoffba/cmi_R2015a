function [image] = CSA_create_maps(V)
%createCSA_maps.m Creates CSA maps from vessel segmentations
% Inputs:
%   V : (3D binary image) vessel segmentation
%
% Outputs: 
%   image : CSA skeleton 

    rng default
    addpath('Z:\CT_Lung\Matlab_Scripts\Emily\BeaumontBatch\Info_Generation\skel2graph3d')
    
    skel = bwskel(V > 0); %Create skeleton of vessels
    mask = bwconncomp(skel, 26); %Break into connected pieces
    
    pieces = mask.PixelIdxList; 

    D = bwdist(V == 0); %distance from each vessel pixel to the background ie radius
    csa = (D.*0.625).^2.*pi;
    
    %Skeleton to fill with CSA 
    image = zeros(size(V));
    
    for obj = 1:length(pieces) %For each connected piece
       %Isolate just this piece on 3D matrix
        tmp = zeros(size(V));
        tmp( pieces{obj}) = true;
        try
            [~,~,links] = Skel2Graph3D(tmp,0);
            for lidx = 1:length(links)
                link = links(lidx).point;
                image(link) = mean(csa(link));
            end
        catch E
        end
    end
    image = single(image);
end