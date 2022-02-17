function [V] = vesselSeg_subj_boxes(ct, lobes)
    lobes = lobes > 0;
    I = double(ct) .*double(lobes);
    I = Normalize(I);
    clearvars ct
    
    V = zeros(size(I));
    
    %Divide I into num_boxes pieces to fit into memory
    [num_boxes, box_size] = calculateBoxSize(I);
    image_size = size(I);
    
    %MTHT3D Parameters
    no = 12; %number of orientations
    beta = 70; 
    alpha = 0.5; 
    c = 15; %parameters for the Vesselness
    alfa = -1/3; %parameters for the Neuritenees
    scale = [0.5,1,1.5,3:2:13];
    
    
    for i = 1:num_boxes
        if i ~= num_boxes
            box = I(:, :, 1+(i-1)*box_size(3):box_size(3)*i);
            [tmp, ~] = MTHT3D(box,scale,no,beta,alpha,c,alfa);
            V(:, :, 1+(i-1)*box_size(3):box_size(3)*i) = tmp;
        else
            box = I(:, :, 1+(i-1)*box_size(3):image_size(3));
             [tmp, ~] = MTHT3D(box,scale,no,beta,alpha,c,alfa);
            V(:, :, 1+(i-1)*box_size(3):image_size(3)) = tmp;
        end
    end
    V = double(V) .* double(lobes);
    
end