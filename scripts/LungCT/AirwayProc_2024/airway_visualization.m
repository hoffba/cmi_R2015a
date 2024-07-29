% airway_visualization.m - basic script to visualize output of large tree

d='R:\CGalban_Lab\LabMembers\AlexanderBell\AirwayProc_2024\data\output\11001N'; % dir w output files
load([d '\Branches.mat'])
load([d '\BranchNodes.mat'])
load([d '\Gen.mat'])

gen2plot=1:10;
plot_big_tree([Branches(:,2:end) Gen],Nodes,4,gen2plot)

load([d '\11001N_tree_data.mat']) % addon to plot lobe surfaces
hold on
for i=1:5
    trisurf(lobe_surfs{i},'FaceColor','b','FaceAlpha',0.05,'EdgeColor','none')
end