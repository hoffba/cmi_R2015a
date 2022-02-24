function [num_boxes, box_size] = calculateBoxSize(I)
    num_boxes = 1;
    box_size = size(I);
    while(prod(box_size) > 140000000)
        box_size(3) = floor(box_size(3)/2);
        num_boxes = num_boxes + 1;
    end
    disp(strcat('number of boxes: ', string(num_boxes)))
end