function plot_tree(ha,B,N,k)
% given branch matrix B ~ [prox,dist,radius] and node matrix [id,x,y,z]
% produce a plot with line thickness set by radii in B
% Update 7/28/20: B appended with generation, Strahler and Horsfield info
% Idea: color lines using one of these features (generation).
% k: 4 ~ generation, 5 ~ Stahler, 6 ~ Horsfield

s = size(B,1);

hold(ha,'on');
grid(ha,'on');
view(ha,[0,1,0]);

maxC=max(B(:,k)); alph=0.9;
for i=1:s
    
    n1=N(ismember(N(:,1),B(i,1)),2:4);
    n2=N(ismember(N(:,1),B(i,2)),2:4);
    
    %plot3([n1(1) n2(1)],[n1(2) n2(2)],[n1(3) n2(3)],'black-',...
    %    'LineWidth',B(i,3))

    colFac=(B(i,k)-1)/(maxC-1);
    colFac=alph*colFac;
    colFac=1-colFac;
    
    plot3(ha,[n1(1) n2(1)],[n1(2) n2(2)],[n1(3) n2(3)],'-',...
        'LineWidth',1,'Color',[1-colFac 0 colFac])
    %B(i,3)
    plot3(ha,n2(1),n2(2),n2(3),'blacko','MarkerFaceColor','black',...
        'MarkerSize',1)
    
end
%plot3(N(:,2),N(:,3),N(:,4),'ro')