function plot_tree(ha,B,Blbl,N,Bind,val_name)
% Inputs:
%   ha          = handle to axes for plotting (if empty a new figure and axes are created)
%   B           = Branches matrix
%   Blbl        = cellstr of labels for columns in B
%   N           = Nodes matrix
%   plotflag    = string describing how to plot
%   Bind        = 1) index for column of B to color the plot
%                 2) values to plot on branches B
%   val_name    = string to title the plot

nB = size(B,1);
cTag = find(strcmp(Blbl,'Tag'));
if ~isempty(cTag)
    nsub = max(B(:,cTag));
end

% % Parse inputs
% p = inputParser();
% addParameter(p,'Lobes',[],@(x)validateattributes(x,'cell',{}));
% addParameter(p,'BranchColor','k',@(x)validateattributes(x,{'numeric'},{'size',[nB,1]}));
% addParameter(p,'Subtree',[],@(x)validateattributes(x,{'numeric'},{'positive','integer','<=',nsub}));
% addParameter(p,'Label','',@(x)validateattributes(x,{'string','char'},{}));




% Prep figure / axes
if isempty(ha)
    hf = figure('Name','Final Visual'); ha = axes(hf);
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