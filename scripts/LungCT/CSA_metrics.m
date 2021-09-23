function [CSA_EXP_A, CSA_EXP_B] = CSA_metrics(image)
%     skel = image > 1.2272;
    skel = image > 0;
    
    desired_edges = zeros(1,11);
    for i = 0:10
        desired_edges(i+1) = 1.2272 + (i)*2;
    end
    
    CSA_EXP_A = NaN;
    CSA_EXP_B = NaN;
    
    [~,~,link] = Skel2Graph3D(skel,0);
    csa_list = zeros(1, length(link));
    for i = 1:length(link)
        point = link(i).point(1);
        csa_list(i) = image(point);
    end
    [vals_per, ~] =histcounts(csa_list, desired_edges, 'Normalization', 'probability');
%     [vals, ~] =histcounts(csa_list, desired_edges);

    centers = zeros(1, length(desired_edges)-1);
    for i = 1:(length(desired_edges)-1)
        centers(i) = mean([desired_edges(i), desired_edges(i+1)]);
    end

    %remove when %of vessels is less than 0.01
    centers = centers(vals_per > 0.01);
%     vals = vals(vals_per > 0.01);
    vals_per = vals_per(vals_per > 0.01);


    try
        m1 = fit(centers', vals_per','exp1');
        CSA_EXP_A = m1.a;
        CSA_EXP_B = m1.b;

    catch E
    end

%         figure; 
%         hold on
%         plot(m2,centers', vals')
%         bar(centers, vals)
%         scatter(centers, vals, 'r.');
%         ylabel('% of Vessels');
%         xlabel('CSA')   
%         title('Distribution of vessel CSA')
%         legend off
%         hold off

end