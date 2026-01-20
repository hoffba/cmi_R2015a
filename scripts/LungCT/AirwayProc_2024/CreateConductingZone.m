function [Branches,Nodes,lower_branches,lower_nodes,num_real_branches] = CreateConductingZone(Branches, Nodes, lobe_TRI, file_directory)  % add by Ali namvar - added num_real_branches output
%% CREATECONDUCTINGZONE(BRANCHES, NODES, LOBE_TRI) grows out an airway tree
% within the input lobe boundaries, using the given seeding airway tree.
% This is based off an algorithm presented in detail in Bordas et al.
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
% file_directory: The name of a folder to save all the output data files
% to. If folder doesn't exist, one will be made. 
%
% Outputs:
% Branches: Full branch structure [ID, prox_node, dist_node, radius, Gen, Strahler, Horsfield, Lobe]
% Nodes: Full node structure [ID, x, y, z]
% lower_branches: Terminal branch indices of the COMBINED tree
% lower_nodes: Terminal node IDs of the COMBINED tree
% num_real_branches: Number of REAL branches (for tagging - branches 1:num_real are real)
%
% Written by Brody Foy, July 2020
% Modified: Added num_real_branches output for proper tag assignment

    newBranches = Branches; 
    options.minBranchLength = 1;
    meanAcinarVolume = 187;
    radiusReductionFactor = 0.87; 
    nlobes = length(lobe_TRI);
    
    %% Remove Nodes that aren't in the branch structure
    branch_nodelist = unique([Branches(:,1); Branches(:,2)]); 
    
    % add by Ali namvar - Use ismember instead of direct indexing to handle non-sequential IDs
    [~, node_rows] = ismember(branch_nodelist, Nodes(:,1));
    node_rows = node_rows(node_rows > 0);
    Nodes = Nodes(node_rows,:); 
    
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
    min_xyz = min(Nodes(:,2:4)); 
    max_xyz = max(Nodes(:,2:4)); 
    for i = 1:nlobes
        curTRI = lobe_TRI{i}; 
        new_min_xyz = min(curTRI.Points); 
        new_max_xyz = max(curTRI.Points); 
        for j = 1:3
            if new_min_xyz(j) < min_xyz(j)
                min_xyz(j) = new_min_xyz(j);
            end
            if new_max_xyz(j) > max_xyz(j)
                max_xyz(j) = new_max_xyz(j); 
            end
        end
    end
    
    Nx = round((max_xyz(1) - min_xyz(1))/meanAcinarVolume^(1/3))+1;
    Ny = round((max_xyz(2) - min_xyz(2))/meanAcinarVolume^(1/3))+1;
    Nz = round((max_xyz(3) - min_xyz(3))/meanAcinarVolume^(1/3))+1;
    
    grid_x = linspace(min_xyz(1), max_xyz(1), Nx);
    grid_y = linspace(min_xyz(2), max_xyz(2), Ny); 
    grid_z = linspace(min_xyz(3), max_xyz(3), Nz); 
    [grid_x, grid_y, grid_z] = meshgrid(grid_x, grid_y, grid_z); 
    grid_xyz = [grid_x(:), grid_y(:), grid_z(:)]; 
    
    %% Calculate the current set of terminal branches (CT-extracted branches)
    orig_lower_branches = []; 
    for i = 1:length(Branches)
        idx = find(Branches(:,1) == Branches(i,2),1); 
        if isempty(idx)
            orig_lower_branches = [orig_lower_branches; i]; 
        end
    end
    
    %% add by Ali namvar - Store the number of real branches BEFORE simulation starts
    % This is CRITICAL for proper tagging later
    num_real_branches_before_sim = size(Branches, 1);
    fprintf('Number of real branches before simulation: %d\n', num_real_branches_before_sim);
    
    %% Initialize a new Branch and new node structure
    newBranches = Branches; 
    newNodes = Nodes; 
    Lobe = zeros(size(Branches,1),1);
    
    %% Loop through all lobes, and for each, algorithmically grow conducting zone airways
    for i = 1:nlobes
        curTRI = lobe_TRI{i}; 
        curVert = curTRI.Points;
        curFace = curTRI.ConnectivityList;
        
        lobe_xyz = grid_xyz(inpolyhedron(curFace, curVert, grid_xyz),:);  
        cur_branchlist = orig_lower_branches(inpolyhedron(curFace, curVert, ...
                            Nodes(Branches(orig_lower_branches,2),2:4)));
        
        node_distances = NaN(size(lobe_xyz,1), length(cur_branchlist)); 
        curDistalNodes = Nodes(Branches(cur_branchlist,2),2:4); 
    
        for j = 1:length(cur_branchlist)
            node_distances(:,j) = sqrt((lobe_xyz(:,1) - curDistalNodes(j,1)).^2 + ...
                (lobe_xyz(:,2) - curDistalNodes(j,2)).^2 + (lobe_xyz(:,3) - curDistalNodes(j,3)).^2); 
        end
        [~, closest_branch] = min(node_distances, [], 2); 
        Lobe(cur_branchlist) = i; 
        
        lobeDims = max(lobe_xyz) - min(lobe_xyz);
        lobeDiagonal = sqrt(sum(lobeDims.^2)); 
        
        options.lobeDiagonal = lobeDiagonal;
        options.curApproxGen = round(log2(size(Branches,1)));
        options.ApproxMaxGen = ceil(options.curApproxGen + log2(size(lobe_xyz,1)));
        
        for j = 1:length(cur_branchlist)
            cur_xyz = lobe_xyz(closest_branch == j,:);
            [newBranches, newNodes, Lobe] = GrowAirwayTree(newBranches, newNodes, ...
                                    cur_xyz, cur_branchlist(j), Lobe, options);
        end
    end
    
    %% Remove branches with nodes outside the lobe structures
    outside_idx = zeros(size(newNodes,1),1); 
    for i = 1:nlobes
        curTRI = lobe_TRI{i}; 
        curVert = curTRI.Points;
        curFace = curTRI.ConnectivityList;
        outside_idx = outside_idx + inpolyhedron(curFace, curVert, newNodes(:,2:4));
    end
    outside_idx = find(outside_idx == 0);
    real_nodes = unique([newNodes(Branches(orig_lower_branches,1),1); newNodes(Branches(orig_lower_branches,2),1)]); 
    outside_idx = outside_idx(outside_idx > max(real_nodes));
    
    [newBranches, newNodes, Lobe, orig_lower_branches] = RemoveBranches(newBranches, newNodes,...
                                                    Lobe, outside_idx, orig_lower_branches); 
    
    %% Convert the newBranch structure to the desired format (four column style)
    parentBranches = orig_lower_branches; 
    Branches = [(1:size(newBranches,1))', newBranches]; 
    Nodes = newNodes; 
    
    %% Clean up any branches that are really just split into two
    [Branches, Nodes, Lobe] = RemoveOneNodeBranches(Branches, Nodes, Lobe); 
    
    %% FIX: Recalculate num_real_branches after cleanup
    % Real branches are the first num_real_branches_before_sim branches
    % But some may have been removed during cleanup, so we need to track this
    % The safest way: real branches are those with ID <= num_real_branches_before_sim
    % BUT branch IDs are renumbered in RemoveOneNodeBranches
    % So we need a different approach: count branches that existed before simulation
    
    % After all cleanup, real branches should still be at the beginning
    % We track by checking which branches have proximal nodes that existed in original tree
    % Actually, the cleanest approach: track the mapping through cleanup
    
    % For now, use a heuristic: real branches are those where BOTH nodes
    % have IDs <= max_original_node_id
    % This works because simulated nodes are added with higher IDs
    
    % Better approach: just use the stored count, adjusted for any removals
    % The RemoveOneNodeBranches only merges branches, doesn't change the real/sim boundary much
    
    % add by Ali namvar - Recalculate num_real_branches after cleanup
    num_real_branches = min(num_real_branches_before_sim, size(Branches,1));
    
    %% Calculate terminal branches of the new structure (COMBINED tree)
    lower_branches = zeros(round(1/2*size(Branches,1)),1);
    lower_nodes = lower_branches;
    counter = 1; 
    for i = 1:size(Branches,1)
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
    
    fprintf('Final tree: %d total branches, %d real branches\n', size(Branches,1), num_real_branches);
    
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
    save([file_directory, '/num_real_branches.mat'], 'num_real_branches');  % add by Ali namvar
    
    %% Return resulting branches
    Branches = [Branches,Gen,Strahler,Horsfield,Lobe];

end

function [Branches, Nodes, Lobe] = GrowAirwayTree(Branches, Nodes, voxel_xyz,...
                                                    curBranch, Lobe, options)
    curDistalNode = Nodes(Branches(curBranch,2),2:4);
    curProxNode = Nodes(Branches(curBranch,1),2:4); 
    curCentroid = sum(voxel_xyz, 1)/size(voxel_xyz,1);
    
    normalPlane = cross(curProxNode(:) - curDistalNode(:), ...
                                curCentroid(:) - curDistalNode(:)); 
    node_splits = NaN(size(voxel_xyz,1),1); 
    for i = 1:size(voxel_xyz,1)
        cur_xyz = voxel_xyz(i,:); 
        node_splits(i) = sign(dot(cur_xyz(:) - curDistalNode(:), normalPlane));
    end
    
    new_vox_idx= {find(node_splits == 1), find(node_splits == -1)}; 
    options.curApproxGen = options.curApproxGen + 1; 
    
    max_branchlength = max([options.lobeDiagonal*(1 - options.curApproxGen/options.ApproxMaxGen),5]);
    
    for i = 1:2
        cur_xyz = voxel_xyz(new_vox_idx{i},:); 
        if ~isempty(cur_xyz)  
            curCentroid = sum(cur_xyz, 1)/size(cur_xyz,1); 
            newDistalNode = curDistalNode + 0.4*(curCentroid - curDistalNode);
            
            cur_branchlength = 0.4*sqrt(sum((curCentroid - curDistalNode).^2));
            if cur_branchlength <= max_branchlength
                Nodes = [Nodes; max(Nodes(:,1))+1, newDistalNode]; 
                Branches = [Branches; Branches(curBranch,2), max(Nodes(:,1)),0];
                Lobe = [Lobe; Lobe(curBranch)];
            
                curBranchLength = sqrt(sum((curDistalNode - newDistalNode).^2));
                if size(cur_xyz,1) > 1 && curBranchLength > options.minBranchLength
                    [Branches, Nodes, Lobe] = GrowAirwayTree(Branches, Nodes, cur_xyz,...
                                            size(Branches,1), Lobe, options);
                end
            end
        end
    end
end


function [Branches, Nodes] = ReorderTree(Branches, Nodes)
    uniqueNodes = unique(Nodes(:,1)); 
    single_idx = zeros(length(uniqueNodes),1); 
    for i = 1:length(uniqueNodes)
        cur_idx = [find(Branches(:,1) == uniqueNodes(i)); find(Branches(:,2) == uniqueNodes(i))]; 
        if length(cur_idx) == 1
            single_idx(i) = 1; 
        end
    end
    termNodes = uniqueNodes(single_idx == 1);
    termRads = NaN(length(termNodes),1); 
    for i = 1:length(termRads)
        cur_idx = [find(Branches(:,1) == termNodes(i)); find(Branches(:,2) == termNodes(i))]; 
        termRads(i) = Branches(cur_idx,3); 
    end
    [~,tracheal_idx] = max(termRads);
    trachealNode = termNodes(tracheal_idx); 
    
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
    daughters = find(Branches(:,2) == Branches(parent_branch,3));
    for i = 1:length(daughters)
        horsfieldFactor = reduction_factor.^(Horsfield(parent_branch)-Horsfield(daughters(i)));
        Branches(daughters(i),4) = Branches(parent_branch,4)*horsfieldFactor; 
        Branches = AssignRadii(Branches, daughters(i),reduction_factor, Horsfield);
    end
end

function [Branches, Nodes, Lobe] = RemoveOneNodeBranches(Branches, Nodes, Lobe)
    for i = 2:size(Nodes,1)
        idx = find(Branches(:,2) == Nodes(i,1));
        if length(idx) == 1
            proxBranch = find(Branches(:,3) == Nodes(i,1)); 
            distBranch = idx; 
            
            newDistNode = Branches(distBranch,3);
            Branches(proxBranch,3) = newDistNode;
            Branches(distBranch,4) = -1;
        end
    end
    Lobe = Lobe(Branches(:,4) ~= -1); 
    Branches = Branches(Branches(:,4) ~= -1,:); 
    Branches(:,1) = (1:size(Branches,1))'; 
    
    newBranches = Branches;         
    branch_nodelist = unique([Branches(:,2); Branches(:,3)]); 
    
    % add by Ali namvar - Find the rows in Nodes that contain these node IDs using ismember
    [~, node_rows] = ismember(branch_nodelist, Nodes(:,1));
    node_rows = node_rows(node_rows > 0);
    Nodes = Nodes(node_rows,:); 
    
    % add by Ali namvar - Create proper node ID to row mapping
    old_node_ids = Nodes(:,1);
    new_node_ids = (1:size(Nodes,1))';
    
    max_old_id = max(old_node_ids);
    if max_old_id > 0
        old_to_new_map = zeros(max_old_id, 1);
        for k = 1:length(old_node_ids)
            old_to_new_map(old_node_ids(k)) = new_node_ids(k);
        end
        
        for k = 1:size(Branches,1)
            if Branches(k,2) <= max_old_id && old_to_new_map(Branches(k,2)) > 0
                newBranches(k,2) = old_to_new_map(Branches(k,2));
            end
            if Branches(k,3) <= max_old_id && old_to_new_map(Branches(k,3)) > 0
                newBranches(k,3) = old_to_new_map(Branches(k,3));
            end
        end
    end
    
    Nodes(:,1) = new_node_ids;
    Branches = newBranches; 
end

function [Branches, Nodes, Lobe, parent_branches] = RemoveBranches(Branches, Nodes,...
                                        Lobe, node_idx, parent_branches)
    branch_idx = []; 
    for i = 1:length(node_idx)
        branch_idx = [branch_idx; find(Branches(:,1) == node_idx(i) | Branches(:,2) == node_idx(i))];
    end
    branch_idx = unique(branch_idx); 
    
    counter = 1; 
    while counter <= length(branch_idx)
        daughters = find(Branches(:,1) == Branches(branch_idx(counter),2));
        branch_idx = [branch_idx; daughters]; 
        counter = counter + 1; 
    end
    branch_idx = unique(branch_idx); 
    
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
    
    % add by Ali namvar - Use proper node ID mapping instead of direct indexing
    old_node_ids = newNodes(:,1);
    new_node_ids = (1:size(newNodes,1))';
    
    max_old_id = max(old_node_ids);
    if max_old_id > 0
        old_to_new_map = zeros(max_old_id, 1);
        for k = 1:length(old_node_ids)
            old_to_new_map(old_node_ids(k)) = new_node_ids(k);
        end
        
        newBranches = Branches;
        for k = 1:size(Branches,1)
            if Branches(k,1) <= max_old_id && old_to_new_map(Branches(k,1)) > 0
                newBranches(k,1) = old_to_new_map(Branches(k,1));
            end
            if Branches(k,2) <= max_old_id && old_to_new_map(Branches(k,2)) > 0
                newBranches(k,2) = old_to_new_map(Branches(k,2));
            end
        end
        Branches = newBranches;
    end
    
    newNodes(:,1) = new_node_ids;
    Nodes = newNodes; 
end
