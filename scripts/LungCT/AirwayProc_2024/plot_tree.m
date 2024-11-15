function plot_tree(ha,B,N,Bval,val_name)
% given branch matrix B ~ [prox,dist,radius] and node matrix [id,x,y,z]
% produce a plot with line thickness set by radii in B
% Update 7/28/20: B appended with generation, Strahler and Horsfield info
% Idea: color lines using one of these features (generation).
% k: 4 ~ generation, 5 ~ Stahler, 6 ~ Horsfield

nB = size(B,1);

% kstr = {'Node1';...
%         'Node2';...
%         'Radius';...
%         'Generation';...
%         'Strahler';...
%         'Horsfield'};

if isempty(ha)
    hf = figure('Name','Final Visual'); ha = axes(hf);
end
hold(ha,'on');
grid(ha,'on');
view(ha,[1,0,0]);
title(ha,val_name);

% Normalize values to [0 1] scale
Bmin = min(min(Bval),0);
Bval = (Bval - Bmin)/(max(Bval) - Bmin);

for i = 1:nB
    
    n = [ N(N(:,1)==B(i,1),2:4);
          N(N(:,1)==B(i,2),2:4) ];
    
    plot3(ha,n(:,1),n(:,2),n(:,3),'-','LineWidth',1,'Color',[Bval(i) 0 1-Bval(i)]);

    plot3(ha,n(2,1),n(2,2),n(2,3),'blacko','MarkerFaceColor','black','MarkerSize',1);
end