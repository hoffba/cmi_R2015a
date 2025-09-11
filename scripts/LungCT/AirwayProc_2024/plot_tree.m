function ha = plot_tree(ha,B,N,C,L_surfs,val_name)
% Inputs:
%   ha          = handle to axes for plotting (if empty a new figure and axes are created)
%   B           = Branches matrix
%   N           = Nodes matrix
%   plotflag    = string describing how to plot
%   Bind        = 1) index for column of B to color the plot
%                 2) values to plot on branches B
%   val_name    = string to title the plot

if nargin<6
    val_name = '';
end
if nargin<5
    L_surfs = [];
end
if nargin<4
    C = [];
end

nB = size(B,1);

% Prep figure / axes
if isempty(ha)
    hf = figure('Name','Airway Tree'); ha = axes(hf);
end
hold(ha,'on');
grid(ha,'on');
view(ha,[1,0,0]);
title(ha,val_name);

% Display lung lobe surfaces
for i = 1:length(L_surfs)
    trisurf(L_surfs{i},'FaceColor','b','FaceAlpha',0.05,'EdgeColor','none');
end

% Normalize values to [0 1] scale
C_flag = length(C)==nB;
if C_flag
    C = (C-min(C)) / (max(C)-min(C));
end

if isempty(C)
    c = [0 0 0];
else
    c = C;
end
for i = 1:nB
    
    n = [ N(N(:,1)==B(i,2),2:4);
          N(N(:,1)==B(i,3),2:4) ];
    
    if C_flag
        c = [C(i) 0 1-C(i)];
    end
    plot3(ha,n(:,1),n(:,2),n(:,3),'-','LineWidth',1,'Color',c);

    plot3(ha,n(2,1),n(2,2),n(2,3),'blacko','MarkerFaceColor','black','MarkerSize',1);
end