function [V] = vesselSeg_subj_boxes(ct, lobes)

    I = single(ct) .* logical(lobes);
    I = Normalize(I);
    I(isnan(I)) = 0;
    
    V = zeros(size(I));
    
    %Divide I into num_boxes pieces to fit into memory
    % image_size = size(I);
    % nblocks = ceil( prod(image_size) / 14e6 );
    % zlim = round(linspace(0,image_size(3),nblocks+1));
%     [num_boxes, box_size] = calculateBoxSize(I);
    
    %MTHT3D Parameters
    no = 12; %number of orientations
    beta = 70; 
    alpha = 0.5; 
    c = 15; %parameters for the Vesselness
    scale = [0.5,1,1.5,3:2:13];

    % Loop over voxels within the lung
    idx = find(lobes);
    for i = numel(idx)
        V(idx) = MTHT3D_BH(I,idx,scale,no,beta,alpha,c);
    end

    
    % for i = 1:nblocks
    %     fprintf('   block %u/%u: %u - %u\n',i,nblocks,zlim(i)+1,zlim(i+1));
    %     zind = (zlim(i)+1):zlim(i+1);
    %     V(:,:,zind) = MTHT3D_BH(I(:,:,zind),scale,no,beta,alpha,c);
    % end
    % V(~lobes) = 0;
end