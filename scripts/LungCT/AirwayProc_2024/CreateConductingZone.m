function CreateConductingZone(Branches, Nodes, lobe_TRI, pxDim, file_directory)
%% CREATECONDUCTINGZONE(BRANCHES, NODES, LOBE_TRI) grows out an airway tree
% within the input lobe boundaries, using the given seeding airway tree.
% This is based off an algorithmic presented in detail in Bordas et al.
% 2015. Development and analysis of patient-based complete conducting
% airways models. 
%
% Inputs:
% Branches: N x 3 array, where each row represents a
% branch, formatted as: [node_1, node_2, radius]
% Nodes: N x 4 array, where each row represents a node: [node_id, x, y, z]
% lobe_TRI: a 5x1 (or 1x5) cell array, where each cell contains a
% triangulation of a lung lobe. Triangulation should be represented as a
% structural array with fields .Points, and .ConnectivityList. 
% pxDim: A 1x3 vector used to convert Node positions to mm units. 
% file_directory: The name of a folder to save all the output data files
% to. If folder doesn't exist, one will be made. 
% Outputs:
% No direct outputs are returned. In the input file_directory folder, the
% following will be saved:
% Branches.mat: An n x 4 array, containing each Branch, structured as:
% [Branch_number, proximal_node, distal_node, radius].
% Nodes.mat: An n x 4 array, containing each Node, structured as:
% [Node_number, x, y, z]
% lower_branches,mat: A list of terminal branch numbers.
% lower_nodes.mat: A list of distal nodes of all terminal branches
% Gen.mat, Horsfield.mat, Strahler.mat: A list of branch generations, Strahler
% orders, and Horsfield orders for all Branches, ordered identically to
% Branches.mat
% Lobe.mat: A list of which lobe each branch is a part of. List will use
% 1..N to denote each lobe, where N is the length of lobe_TRI. 0 will be
% used to denote any branches not within a lobe
%
% Written by Brody Foy, July 2020

newBranches = Branches; 
options.minBranchLength = 1;
meanAcinarVolume = 187; %This will determine (roughly) how many branches the output tree has.
radiusReductionFactor = 0.87; 

%% Rescale Nodes structure 
for i = 1:3
    Nodes(:,i+1) = Nodes(:,i+1)*pxDim(i); 
end

%% Need to set max branch length, using T = max(Vb - Dl*n, 5mm)
% where Vb = size of the diagonal of the bounding box of the lobe
% Dl = Vb/N, where N is max num of gens, and n is the current gen

%% Remove Nodes that aren't in the branch structure
branch_nodelist = unique([Branches(:,1); Branches(:,2)]); 
Nodes = Nodes(branch_nodelist,:); 

%% First re-align Nodes (so that there aren't missing nodes)
for i = 1:size(Nodes,1)
    new_idx = size(Nodes,1) - i + 1;
    cur_idx = Nodes(end-i+1,1); 
    Nodes(end-i+1,1) = new_idx; 
    newBranches(Branches(:,1) == cur_idx,1) = new_idx;
    newBranches(Branches(:,2) == cur_idx,2) = new_idx; 
end
Branches = newBranches; 

%% Next reorder the branches structure (so that first node is always proximal node)
[Branches, Nodes] = ReorderTree(Branches, Nodes);

N1 = 1; 
N2 = Branches(1,1); 

Nodes(Nodes(:,1) == N1,1) = -1; 
Nodes(Nodes(:,1) == N2,1) = N1; 
Nodes(Nodes(:,1) == -1,1) = N2; 

Branches(Branches(:,1) == N1,1) = -1; 
Branches(Branches(:,1) == N2,1) = N1; 
Branches(Branches(:,1) == -1,1) = N2; 

Branches(Branches(:,2) == N1,2) = -1; 
Branches(Branches(:,2) == N2,2) = N1; 
Branches(Branches(:,2) == -1,2) = N2; 

[~, sort_idx] = sort(Nodes(:,1)); 
Nodes = Nodes(sort_idx,:); 


%% Create an even grid of xyz points, to grow airway tree towards

% First calculate dimensions of lung (from lobe triangulations)
min_xyz = min(Nodes(:,2:4)); 
max_xyz = max(Nodes(:,2:4)); 
for i = 1:5
    curTRI = lobe_TRI{i}; 
    new_min_xyz = min(curTRI.Points).*pxDim; 
    new_max_xyz = max(curTRI.Points).*pxDim; 
    for j = 1:3
        if new_min_xyz(j) < min_xyz(j)
            min_xyz(j) = new_min_xyz(j);
        end
        if new_max_xyz(j) > max_xyz(j)
            max_xyz(j) = new_max_xyz(j); 
        end
    end
end

% Create a 3D grid of points, defining voxel corners, such that each voxel
% has volume: meanAcinarVolume (in mm^3)
Nx = round((max_xyz(1) - min_xyz(1))/meanAcinarVolume^(1/3))+1;
Ny = round((max_xyz(2) - min_xyz(2))/meanAcinarVolume^(1/3))+1;
Nz = round((max_xyz(3) - min_xyz(3))/meanAcinarVolume^(1/3))+1;

grid_x = linspace(min_xyz(1), max_xyz(1), Nx);
grid_y = linspace(min_xyz(2), max_xyz(2), Ny); 
grid_z = linspace(min_xyz(3), max_xyz(3), Nz); 
[grid_x, grid_y, grid_z] = meshgrid(grid_x, grid_y, grid_z); 
grid_xyz = [grid_x(:), grid_y(:), grid_z(:)]; 

%% Calculate the current set of terminal branches (CT-extracted branches)
lower_branches = []; 
for i = 1:length(Branches)
    idx = find(Branches(:,1) == Branches(i,2),1); 
    if isempty(idx)
        lower_branches = [lower_branches; i]; 
    end
end

%% Initialize a new Branch and new node structure
newBranches = Branches; 
newNodes = Nodes; 
Lobe = zeros(size(Branches,1),1); %Holds which Lobe a branch is in

%% Loop through the 5 lobes, and for each, algorithmically grow conducting 
% zone airways
for i = 1:length(lobe_TRI)
    curTRI = lobe_TRI{i}; 
    curVert = curTRI.Points;
    %% Rescale current lobe dimensions, into mm
    for j = 1:3
        curVert(:,j) = curVert(:,j)*pxDim(j); 
    end
    
    curFace = curTRI.ConnectivityList;
    % Calculate which points in the overall grid are inside the given lobe
    lobe_xyz = grid_xyz(inpolyhedron(curFace, curVert, grid_xyz),:);  
    % Terminal branches that are in the current lobe
    cur_branchlist = lower_branches(inpolyhedron(curFace, curVert, ...
                        Nodes(Branches(lower_branches,2),2:4)));
    
    %% Assign each node (in cur_xyz) to its nearest branch
    node_distances = NaN(size(lobe_xyz,1), length(cur_branchlist)); 
    curDistalNodes = Nodes(Branches(cur_branchlist,2),2:4); 

    for j = 1:length(cur_branchlist)
        node_distances(:,j) = sqrt((lobe_xyz(:,1) - curDistalNodes(j,1)).^2 + ...
            (lobe_xyz(:,2) - curDistalNodes(j,2)).^2 + (lobe_xyz(:,3) - curDistalNodes(j,3)).^2); 
    end
    [~, closest_branch] = min(node_distances, [], 2); 
    Lobe(cur_branchlist) = i; 
    
    %% Calculate lobe bounding volume
    lobeDims = max(lobe_xyz) - min(lobe_xyz);
    lobeDiagonal = sqrt(sum(lobeDims.^2)); 
    
    options.lobeDiagonal = lobeDiagonal;
    options.curApproxGen = round(log2(size(Branches,1)));
    options.ApproxMaxGen = ceil(options.curApproxGen + log2(size(lobe_xyz,1)));
    
    
    %% Now call recurring function to grow airways out from each of the
    % terminal branches, to the set of associated points
    for j = 1:length(cur_branchlist)
        cur_xyz = lobe_xyz(closest_branch == j,:);
        [newBranches, newNodes, Lobe] = GrowAirwayTree(newBranches, newNodes, ...
                                cur_xyz, cur_branchlist(j), Lobe, options);
    end
end

%% Remove branches with nodes outside the lobe structures
outside_idx = zeros(size(newNodes,1),1); 
for i = 1:length(lobe_TRI)
    curTRI = lobe_TRI{i}; 
    curVert = curTRI.Points;
    %% Rescale current lobe dimensions, into mm
    for j = 1:3
        curVert(:,j) = curVert(:,j)*pxDim(j); 
    end
    curFace = curTRI.ConnectivityList;
    outside_idx = outside_idx + inpolyhedron(curFace, curVert, newNodes(:,2:4));
end
outside_idx = find(outside_idx == 0); % Not in one of the Lobes
% Need to compare to CT extracted branches (as many of those will be
% outside Lobes, due to being too proximal for the lobe. 
real_nodes = unique([newNodes(Branches(lower_branches,1),1); newNodes(Branches(lower_branches,2),1)]); 
outside_idx = outside_idx(outside_idx > max(real_nodes));

[newBranches, newNodes, Lobe, lower_branches] = RemoveBranches(newBranches, newNodes,...
                                                Lobe, outside_idx, lower_branches); 

%% Convert the newBranch structure to the desired format (four column style)
parentBranches = lower_branches; 
Branches = [(1:size(newBranches,1))', newBranches]; 
Nodes = newNodes; 

%% Clean up any branches that are really just split into two
[Branches, Nodes, Lobe] = RemoveOneNodeBranches(Branches, Nodes, Lobe); 

%% Calculate terminal branches of the new structure
lower_branches = zeros(round(1/2*size(Branches,1)),1);
lower_nodes = lower_branches;
counter = 1; 
for i = 1:length(Branches)
    idx = find(Branches(:,2) == Branches(i,3), 1, 'first');
    if isempty(idx)
        lower_branches(counter) = i; 
        lower_nodes(counter) = Branches(i,3);
        counter = counter + 1; 
    end
end

lower_branches = lower_branches(lower_branches ~= 0);
lower_nodes = lower_nodes(lower_nodes ~= 0); 

%% Calculate Generation, Horsfield and Strahler numbers of the structure
[Gen, Horsfield, Strahler] = Get_Gen_Horsfield(Branches, lower_branches);

%% Now reassign airway radii (for the algorithmic airways) using Bordas et al formula
for i = 1:length(parentBranches)
    daughters = find(Branches(:,2) == Branches(parentBranches(i),3));
    if ~isempty(daughters)
        horsfieldFactor = radiusReductionFactor.^(Horsfield(parentBranches(i))-Horsfield(daughters));
        Branches(daughters,4) = Branches(parentBranches(i),4)*horsfieldFactor; 
        for j = 1:length(daughters)
            Branches = AssignRadii(Branches, daughters(j), radiusReductionFactor, Horsfield); 
        end
    end
end


%% Save output files to the given directory
if ~exist(file_directory, 'dir')
    mkdir(file_directory);
end
save([file_directory, '/Branches.mat'], 'Branches');
save([file_directory, '/BranchNodes.mat'], 'Nodes');
save([file_directory, '/lower_nodes.mat'], 'lower_nodes');
save([file_directory, '/lower_branches.mat'], 'lower_branches');
save([file_directory, '/Strahler.mat'], 'Strahler');
save([file_directory, '/Gen.mat'], 'Gen');
save([file_directory, '/Horsfield.mat'], 'Horsfield');
save([file_directory, '/Lobe.mat'], 'Lobe'); 



end

function [Branches, Nodes, Lobe] = GrowAirwayTree(Branches, Nodes, voxel_xyz,...
                                                    curBranch, Lobe, options)
%% [BRANCHES, NODES, LOBE] = GROWAIRWAYTREE(BRANCHES, NODES, VOXEL_XYZ, CURBRANCH,
% LOBE, OPTIONS) grows two new branches out from the curBranch, using the 
% seed voxels (voxel_xyz), starting from the current Branch and Node
% structure. Branches are grown by splitting the voxel grid into 2 (based
% on the plane from the branch to the voxel centroid). Lobe is a structure
% denoting which Lobe all branches are in. Algorithm is based off of paper
% by Bordas et al. 2015, PLoS One. ('Development and analysis of
% patient-based...')
%
% Written by Brody Foy, July 2020


%% Calculate centroid of voxels
curDistalNode = Nodes(Branches(curBranch,2),2:4);
curProxNode = Nodes(Branches(curBranch,1),2:4); 
curCentroid = sum(voxel_xyz, 1)/size(voxel_xyz,1);

%% Calculate the normal to the plane defined by the distal, proximal and 
% centroid nodes (using a cross product)
normalPlane = cross(curProxNode(:) - curDistalNode(:), ...
                            curCentroid(:) - curDistalNode(:)); 
% Calculate which side of the plane each point (in voxel_xyz) is on
node_splits = NaN(size(voxel_xyz,1),1); 
for i = 1:size(voxel_xyz,1)
    % The dot product of the vector from each point to the plane, against
    % the plane normal vector will have sign 1 or -1 for opposing sides of
    % the plane
    cur_xyz = voxel_xyz(i,:); 
    node_splits(i) = sign(dot(cur_xyz(:) - curDistalNode(:), normalPlane));
end

new_vox_idx= {find(node_splits == 1), find(node_splits == -1)}; 
options.curApproxGen = options.curApproxGen + 1; 

% Max branch length is from formula in Bordas et al. 2015
max_branchlength = max([options.lobeDiagonal*(1 - options.curApproxGen/options.ApproxMaxGen),5]);

%% Loop through each subset of voxels
for i = 1:2
    cur_xyz = voxel_xyz(new_vox_idx{i},:); 
    if ~isempty(cur_xyz)  
        curCentroid = sum(cur_xyz, 1)/size(cur_xyz,1); 
        newDistalNode = curDistalNode + 0.4*(curCentroid - curDistalNode);
        
        %% Compare to max branch length, only add branch if below
        cur_branchlength = 0.4*sqrt(sum((curCentroid - curDistalNode).^2));
        if cur_branchlength <= max_branchlength
            % Add branch, and new distal node (and assign lobe number to
            % branch)
            Nodes = [Nodes; max(Nodes(:,1))+1, newDistalNode]; 
            Branches = [Branches; Branches(curBranch,2), max(Nodes(:,1)),0];
            Lobe = [Lobe; Lobe(curBranch)];
        
            curBranchLength = sqrt(sum((curDistalNode - newDistalNode).^2));
            % If there is more than 1 voxel connected to the current
            % branch, and the branch isn't too small, continue to grow.
            if size(cur_xyz,1) > 1 && curBranchLength > options.minBranchLength
                [Branches, Nodes, Lobe] = GrowAirwayTree(Branches, Nodes, cur_xyz,...
                                        size(Branches,1), Lobe, options);
            end
        end
    end
end

end


function [Branches, Nodes] = ReorderTree(Branches, Nodes)
% [BRANCHES, NODES] = REORDER_TREE(BRANCHES, NODES) reorders the input tree
% structure (defined by the branch structure, and node structural, such
% that the output branch structure is ordered: [Proximal node, distal node,
% branch radius], and the first row is always the trachea
%
% Branch structure is ordered: [node_1, node_2, radius]
% Node structure is ordered: [node_idx, x, y, z]
%
% Written by Brody Foy, July 2020

%% First identify the tracheal node (only appears once in branch structure, 
% and has a large radius). 
uniqueNodes = unique(Nodes(:,1)); 
single_idx = zeros(length(uniqueNodes),1); 
% First find all nodes only appearing once in branch structure
for i = 1:length(uniqueNodes)
    cur_idx = [find(Branches(:,1) == uniqueNodes(i)); find(Branches(:,2) == uniqueNodes(i))]; 
    if length(cur_idx) == 1
        single_idx(i) = 1; 
    end
end
termNodes = uniqueNodes(single_idx == 1);
% Calculate radii of all terminal nodes
termRads = NaN(length(termNodes),1); 
for i = 1:length(termRads)
    cur_idx = [find(Branches(:,1) == termNodes(i)); find(Branches(:,2) == termNodes(i))]; 
    termRads(i) = Branches(cur_idx,3); 
end
[~,tracheal_idx] = max(termRads); %Trachea has the largest radii of a terminal node
trachealNode = termNodes(tracheal_idx); 


%%  Then reorder the tracheal branch, and recursively call the realignment function
newBranches = zeros(size(Branches,1),3); 
cur_idx = find(Branches(:,1) == trachealNode);
if ~isempty(cur_idx)
    distalNode = Branches(cur_idx,2); 
else
    distalNode = Branches(Branches(:,2) == trachealNode,1); 
end

newBranches(1,:) = [trachealNode, distalNode, termRads(tracheal_idx)]; 
Branches = AlignBranches(newBranches, Branches, trachealNode, distalNode);


end

function newBranches = AlignBranches(newBranches, Branches, proxNode, distalNode)
%% NEWBRANCHES = ALIGNBRANCHES(NEWBRANCHES, BRANCHES, PROXNODE, DISTALNODE) 
% realigns the branch structure, such that each row of the Branch array is
% ordered: [prox_node, distal_node, radius]. The original structure is
% ordered without specific node direction (i.e. proximal node can be in
% either column 1 or 2). Function operates recursively, starting from the
% trachea. 
%
% Written by Brody Foy, July 2020
newBranch_idx = [find(Branches(:,1) == distalNode); find(Branches(:,2) == distalNode)];

for i = 1:length(newBranch_idx)
   counter = find(newBranches(:,1) == 0,1); 
    
   newRad = Branches(newBranch_idx(i), 3); 
   newDistalNode = Branches(newBranch_idx(i),1:2); 
   newDistalNode = newDistalNode(newDistalNode ~= distalNode & newDistalNode ~= proxNode); 
   if ~isempty(newDistalNode)
       newBranches(counter,:) = [distalNode, newDistalNode, newRad]; 
       newBranches = AlignBranches(newBranches, Branches, distalNode, newDistalNode); 
   end
    
end

end

function Branches = AssignRadii(Branches, parent_branch, reduction_factor, Horsfield)
%% BRANCHES = ASSIGNRADII(BRANCHES, PARENT_BRANCH, REDUCTION_FACTOR, HORSFIELD) 
% assigns radii to the Branches in the given Branch structure, using a recursive
% algorithm, where daughters have a radius that is a constant multiple of
% their parent branch radius, to the power of difference between horsfield numbers. 
% 
% Written by Brody Foy, July 2020
daughters = find(Branches(:,2) == Branches(parent_branch,3));
for i = 1:length(daughters)
    horsfieldFactor = reduction_factor.^(Horsfield(parent_branch)-Horsfield(daughters(i)));
    Branches(daughters(i),4) = Branches(parent_branch,4)*horsfieldFactor; 
    Branches = AssignRadii(Branches, daughters(i),reduction_factor, Horsfield);
end

end

function [Branches, Nodes, Lobe] = RemoveOneNodeBranches(Branches, Nodes, Lobe)
%% [BRANCHES, NODES, LOBE] = REMOVEONENODEBRANCHES(BRANCHES, NODES, LOBE) cleans up the
% branch structure by replacing any branches that are split into two, by a
% single branch (with averaged radii). 

for i = 2:length(Nodes) % skip tracheal node
    idx = find(Branches(:,2) == Nodes(i,1)); %Find all branches associated with that proximal node
    if length(idx) == 1
        proxBranch = find(Branches(:,3) == Nodes(i,1)); 
        distBranch = idx; 
        
        newDistNode = Branches(distBranch,3);
        Branches(proxBranch,3) = newDistNode;
        Branches(distBranch,4) = -1; %Mark radius as -1, to flag for removal (at end)
    end
end
Lobe = Lobe(Branches(:,4) ~= -1); 
Branches = Branches(Branches(:,4) ~= -1,:); 
Branches(:,1) = [1:size(Branches,1)]'; 

%% Now clean up nodes (remove now defunt nodes, and reorder)
newBranches = Branches;         
branch_nodelist = unique([Branches(:,2); Branches(:,3)]); 
Nodes = Nodes(branch_nodelist,:); 
% Re-align Nodes (so that there aren't missing nodes)
for i = 1:size(Nodes,1)
    new_idx = size(Nodes,1) - i + 1;
    cur_idx = Nodes(end-i+1,1); 
    Nodes(end-i+1,1) = new_idx; 
    newBranches(Branches(:,2) == cur_idx,2) = new_idx;
    newBranches(Branches(:,3) == cur_idx,3) = new_idx; 
end
Branches = newBranches; 

end

function [Branches, Nodes, Lobe, parent_branches] = RemoveBranches(Branches, Nodes,...
                                        Lobe, node_idx, parent_branches)
%% [BRANCHES, NODES, LOBE, PARENT_BRANCHES] = REMOVEBRANCHES(BRANCHES, NODES, 
% LOBE, NODE_IDX, PARENT_BRANCHES) removes the nodes indicated by Node_idx, 
% as well as any branches containing those nodes, and any daughters of those
% branches. It then reorders the branch and node structures to account for this. 
% parent_branches stores the branch numbers of the original CT airways.
% This is adjusted as branches are removed to reflect the new numbering.
%
% Written by Brody H Foy, Oct 2020

%% First find all affected branches (either with proximal or distal nodes in node_idx)
branch_idx = []; 
for i = 1:length(node_idx)
    branch_idx = [branch_idx; find(Branches(:,1) == node_idx(i) | Branches(:,2) == node_idx(i))];
end
branch_idx = unique(branch_idx); 

%% Now find all daughters of those branches 
counter = 1; 
while counter <= length(branch_idx)
    daughters = find(Branches(:,1) == Branches(branch_idx(counter),2));
    branch_idx = [branch_idx; daughters]; 
    counter = counter + 1; 
end
branch_idx = unique(branch_idx); 

%% Remove all Branches (and associated Lobes), and Nodes, matching to those
node_idx = unique([Branches(branch_idx,2)]); 
[~, new_node_idx] = setdiff(Nodes(:,1), node_idx); 
newNodes = Nodes(new_node_idx,:); 
[~, new_branch_idx] = setdiff(1:size(Branches,1), branch_idx);

newBranches = Branches(new_branch_idx,:); 
for i = 1:length(parent_branches)
    newBranch = find(newBranches(:,1) == Branches(parent_branches(i),1) &...
                     newBranches(:,2) == Branches(parent_branches(i),2));
    parent_branches(i) = newBranch; 
end
Branches = newBranches;
Lobe = Lobe(new_branch_idx); 

%% Re-number all nodes and branches, to not have missing nodes, etc. 
for i = size(newNodes,1):-1:1
    old_idx = newNodes(i,1);
    new_idx = i; 
    
    newBranches(Branches(:,1) == old_idx,1) = new_idx; 
    newBranches(Branches(:,2) == old_idx,2) = new_idx; 
end
newNodes(:,1) = [1:size(newNodes,1)]'; 
Branches = newBranches; 
Nodes = newNodes; 
end
