function [vol, numVessels] = CSA_range_metrics(V_idx, image)
    vol = 0;
    numVessels = 0;
    [~,~,link] = Skel2Graph3D(V_idx,0);
    numVessels = numVessels + length(link);
    for l = 1:length(link)
        try
            vol = vol + image(link(l).point(1))*length(link(l).point)*0.625^3;
        catch E
        end
    end
end
%imreconstruct retrieves all connected pieces instead of just the skeleton