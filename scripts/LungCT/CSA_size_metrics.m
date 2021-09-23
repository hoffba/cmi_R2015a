function [NUM_VESSELS, NUM_COMPONENTS, NUM_ENDPOINTS, AVG_ENDPOINTS_PER_COMPONENT] = CSA_size_metrics(image)
    
%     skel = image > 1.2272;
    skel = image > 0;
    [~,node,link] = Skel2Graph3D(skel,0);
    tbl = struct2table(node);
    NUM_ENDPOINTS = sum(tbl.ep);
    NUM_VESSELS = length(link);

    mask = bwconncomp(skel, 26);
    NUM_COMPONENTS  = length(mask.PixelIdxList);
    AVG_ENDPOINTS_PER_COMPONENT = NUM_ENDPOINTS/NUM_VESSELS - 2;
end