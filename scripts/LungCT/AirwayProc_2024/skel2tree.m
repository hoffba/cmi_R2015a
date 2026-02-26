function [B,N,points] = skel2tree(CL,voxsz)


% Generate rough tree with nodes and links
[~,node,link] = Skel2Graph3D(CL,0);

% Re-format data into the desirerd Branch (B) and Node (N) matrices
B = [link.n1 ; link.n2]';
N = [(1:numel(node))',[node.comx;node.comy;node.comz]' .* voxsz];

% Extract centerline points assigned to each branch
nlink = numel(link);
points = cell(nlink,1);
for i = 1:nlink
    [xx,yy,zz] = ind2sub(size(CL),link(i).point');
    points{i} = [xx,yy,zz] .* voxsz;
end