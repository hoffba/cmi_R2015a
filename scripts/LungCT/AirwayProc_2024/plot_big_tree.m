function plot_big_tree(B,N,k,k2plot)
% given branch matrix B ~ [prox,dist,radius] and node matrix [id,x,y,z]
% produce a plot with line thickness set by radii in B
% Update 7/28/20: B appended with generation, Strahler and Horsfield info
% Idea: color lines using one of these features (generation).
% k: 4 ~ generation, 5 ~ Stahler, 6 ~ Horsfield
% k2plot ~ array of values of feature k to be accepted for plotting

s=size(B,1);

hold on
grid on

%maxC=max(B(:,k));
maxC=max(k2plot); 
alph=0.9;

for i=1:s
    
    if ~ismember(B(i,k),k2plot)
        continue
    end
    
    n1=N(ismember(N(:,1),B(i,1)),2:4);
    n2=N(ismember(N(:,1),B(i,2)),2:4);

    colFac=(B(i,k)-1)/(maxC-1);
    colFac=alph*colFac;
    colFac=1-colFac;
    
    plot3([n1(1) n2(1)],[n1(2) n2(2)],[n1(3) n2(3)],'-',...
        'LineWidth',2*B(i,3),'Color',[1-colFac 0 colFac])

end