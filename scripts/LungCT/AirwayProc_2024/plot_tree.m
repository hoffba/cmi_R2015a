function plot_tree(ha,B,N,Bind,val_name)
% Inputs:
%   ha = handle to axes for plotting (if empty a new figure and axes are created)
%   B = Branches matrix
%   N = Nodes matrix
%   Bind = 1) index for column of B to color the plot
%          2) values to plot on branches B
%   val_name = string to title the plot

nB = size(B,1);

if isempty(ha)
    hf = figure('Name','Final Visual'); ha = axes(hf);
end
hold(ha,'on');
grid(ha,'on');
view(ha,[1,0,0]);
title(ha,val_name);

% Normalize values to [0 1] scale
if numel(Bind)==1
    Bmin = min(min(B(:,Bind)),0);
    Bval = (B(:,Bind) - Bmin)/(max(B(:,Bind)) - Bmin);
else
end

for i = 1:nB
    
    n = [ N(N(:,1)==B(i,1),2:4);
          N(N(:,1)==B(i,2),2:4) ];
    
    plot3(ha,n(:,1),n(:,2),n(:,3),'-','LineWidth',1,'Color',[Bval(i) 0 1-Bval(i)]);

    plot3(ha,n(2,1),n(2,2),n(2,3),'blacko','MarkerFaceColor','black','MarkerSize',1);
end