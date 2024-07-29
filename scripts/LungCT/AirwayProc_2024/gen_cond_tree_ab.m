% gen_cond_tree_ab.m
% using tree prep data, implement Brody's code to generate cond tree struct
pID='10123W';
inD='C:\Users\alejbell\Desktop\Tree Growth Pre-data work\data\output\old\';
outD='C:\Users\alejbell\Desktop\Tree Growth Pre-data work\testing';

load([inD pID '\' pID '_tree_data.mat'])

CreateConductingZone(B(:,1:3), N, lobe_surfs,[1 1 1],outD)

cd(outD);

load('Branches.mat')
load('BranchNodes.mat')
load('Strahler.mat')
load('Gen.mat')
load('Horsfield.mat')
load('Lobe.mat')
load('lower_branches.mat')
load('lower_nodes.mat')

cd('C:\Users\alejbell\Desktop\Tree Growth Pre-data work\workspace')

B=[Branches(:,2:end) Gen Strahler Horsfield];
N=Nodes;

plot_big_tree(B,N,4,1:12)
hold on
for i=1:5
    trisurf(lobe_surfs{i},'FaceColor','b','FaceAlpha',0.05,'EdgeColor','none')
end
axis equal