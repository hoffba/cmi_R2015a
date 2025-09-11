function iowa_preprocess(datapath,savepath)
% Takes airways data from Iowa and saves it in the format required by MiTAP

ID = dir(datapath);
ID(1:2) = [];
ID(~[ID.isdir]) = [];
ID = {ID.name}';
N = numel(ID);

try

for i = 1:N

    fprintf('Preprocessing %s ... \n',ID{i})

    procdir = fullfile(savepath,[ID{i},'_']);
    if ~isfolder(procdir)
        mkdir(procdir);
    end

    airdir = fullfile(procdir,[ID{i},'_.AirwayProc']);
    if ~isfolder(airdir)
        mkdir(airdir);
    end

    % Check that all necessary files are there
    fn_B = fullfile(datapath,ID{i},[ID{i},'_Branches.csv']);
    fn_N = fullfile(datapath,ID{i},[ID{i},'_Junctions.csv']);
    fn_airways = fullfile(datapath,ID{i},'Airway_Mask.nii.gz');
    fn_CL = fullfile(datapath,ID{i},'Airway_Centerline.nii.gz');

    % File to save branches and nodes to
    fn_realtree = fullfile(airdir,'RealTree.mat');

    if isfile(fn_B) && isfile(fn_N) && isfile(fn_CL) && isfile(fn_airways)

        fprintf('    Compiling branches and nodes ...\n');

        [A,~,fov,~,~] = readNIFTI(fn_airways);
        d = size(A);
        % voxsz = fov./d;

        % Branches:
        B = table2array(readtable(fn_B));
        B(:,1) = []; % Don't need a branch ID
        nB = size(B,1);

        % Nodes - need to convert to MiTAP xyz coordinates
        N = readtable(fn_N);
        N = [N.Node,N.X+1,N.Y+1,N.Z+1];
        N_idx = sub2ind(d,N(:,2),N(:,3),N(:,4));
        % nN = size(N,1);

        % Prep skeleton for radius calculations
        % [~,A_D_idx] = bwdist(~A,'euclidean'); % to use for radii approximation
        skel = logical(readNIFTI(fn_CL));
        [~, skel_D_idx] = bwdist(skel);
        % skel_linear_indices = find(skel);
        % num_skel_points = length(skel_linear_indices);
        % linear_to_node_id = zeros(d);
        % linear_to_node_id(skel_linear_indices) = 1:num_skel_points;

        % Quick graph of skeleton
        [~,node,link] = Skel2Graph3D(skel,0);
        link_point = {link.point}';
        link_point_all_idx = [link_point{:}]';
        [rr,cc,ss] = ind2sub(d,link_point_all_idx);
        link_point_all = [rr,cc,ss];

        % Make sure all nodes are on skeleton graph
        missing_point = find(~ismember(N_idx,link_point_all_idx));
        k = dsearchn(link_point_all,N(missing_point,2:4));
        for ik = 1:length(k)

            % Find link associated with nearest point
            k_idx = link_point_all_idx(k(ik));
            klink = find(cellfun(@(x)ismember(k_idx,x),link_point),1);

            % Find location relative to link points list
            point_list_idx = link_point{klink}';
            np = length(point_list_idx);
            [rr,cc,ss] = ind2sub(d,point_list_idx);
            point_list = [rr,cc,ss];
            [~,dist] = dsearchn(N(missing_point(ik),2:4),point_list);
            pdiff = sqrt(sum(diff(point_list,1).^2,2));
            kpt = find(dist==min(dist));
            mp_idx = N_idx(missing_point(ik));
            % Decide where to put it
            if kpt == 1
                if dist(2)<pdiff(1)
                    btwn = [1 2];
                else
                    btwn = [0 1];
                end
            elseif kpt == numel(dist)
                if dist(end-1)<pdiff(end)
                    btwn = [np-1 np];
                else
                    btwn = [np np+1];
                end
            else
                kmin = dist==min(dist);
                knext = find(dist==min(dist(~kmin)),1);
                kmin = find(kmin,1);
                btwn = [min(kmin,knext),max(kmin,knext)];
            end
            link_point{klink} = [ point_list_idx(1:btwn(1)) ; mp_idx ; point_list_idx(btwn(2):end) ]';

        end
    
        % Find voxels just outside airways for radius calculation
        A_border_idx = find(imdilate(A,strel('sphere',1)) & ~A);
        [rr,cc,ss] = ind2sub(d,A_border_idx);
        A_border_point = [rr,cc,ss];

        % Loop over branches and calculate radii
        fprintf('    Calculating branch mean radii ...\n')
        for iB = 1:nB

            % start node
            ind = N(:,1)==B(iB,1);
            start_node = N(ind,2:4);
            start_node_idx = N_idx(ind);
            % Check if node is on the skeleton
            if ~skel(start_node(1), start_node(2), start_node(3))
                idx = sub2ind(d, start_node(1), start_node(2), start_node(3));
                start_node_idx = skel_D_idx(idx);
                % [start_node(1), start_node(2), start_node(3)] = ind2sub(d, start_node_idx);
            end

            % end node
            ind = N(:,1)==B(iB,2);
            end_node =   N(ind,2:4);
            end_node_idx = N_idx(ind);
            % Check if node is on the skeleton
            if ~skel(end_node(1), end_node(2), end_node(3))
                idx = sub2ind(d, end_node(1), end_node(2), end_node(3));
                end_node_idx = skel_D_idx(idx);
            end

            % Find SkelGraph links containing start and end nodes
            ind = cellfun(@(x)any(ismember([start_node_idx,end_node_idx],x)),link_point);

            % Create graph with link points
            gpoint = unique([link_point{ind}]);
            ngp = length(gpoint);
            G = graph;
            G = addnode(G,ngp);
            for igp = 1:ngp
                [r, c, p] = ind2sub(d, gpoint(igp));
    
                % Use 26-connectivity for neighbors.
                [neighbor_r, neighbor_c, neighbor_p] = ndgrid(r-1:r+1, c-1:c+1, p-1:p+1);
                neighbors_idx = sub2ind(d, neighbor_r(:), neighbor_c(:), neighbor_p(:));
                neighbors_idx(14) = [];  % Exclude point itself

                % Find neighbors in set of skeleton points
                ind = find(ismember(gpoint,neighbors_idx));
    
                for i_neighbor = ind'
                    % Add edge only once to avoid duplicates.
                    G = addedge(G, igp, i_neighbor);
                end
            end

            % Find the shortest path between start and end nodes
            % Find the graph node IDs for the start and end points.
            start_id = find(start_node_idx==gpoint);
            end_id = find(end_node_idx==gpoint);
            path_nodes = shortestpath(G, start_id, end_id);
            
            % Convert path back to voxel coordinates
            path_linear_indices = gpoint(path_nodes);
            [path_rows, path_cols, path_planes] = ind2sub(d, path_linear_indices');
            path_points = [path_rows, path_cols, path_planes];
    
            % Calculate mean radius of branch
            % - only calculate middle section of branch
            if ~isempty(path_points)

                % Use points with margin away from nodes
                np = size(path_points,1);
                c = (np+1)/2;
                lim = round([c-np/4 c+np/4]);
                path_points = path_points(lim(1):lim(2),:);
                % np = size(path_points,1);

                % Find distance to nearest airway border
                [~,dist] = dsearchn(A_border_point,path_points);

                B(iB,3) = mean(dist); % approximate radius in mm from mean of sampled midpoints
            else
                B(iB,3) = nan;
            end
        end

        % Save tree to file for MiTAP to load
        B_label = {'N_Prox','N_Dist','Radius'};
        save(fn_realtree,'N','B','B_label');
    end
    
    % Copy certain files to the processing directory
    fn_copy = {'CT_Expiratory',         'exp';...
               'CT_Inspiratory',        'ins';...
               'Lung_Lobes_Expiratory', 'exp.label';...
               'Lung_Lobes_Inspiratory','ins.label';...
               'Airway_Mask',           'airways'};
    for j = 1:size(fn_copy,1)
        fprintf('    Copying file: %s to %s\n',fn_copy{j,:});
        fn = fullfile(datapath,ID{i},[fn_copy{j,1},'.nii.gz']);
        if isfile(fn)
            copyfile(fn,fullfile(procdir,[ID{i},'_.',fn_copy{j,2},'.nii.gz']),CopyLinkBehavior="resolve");
        end
    end

    fprintf('  done\n')

end

catch err
    disp('breakpoint')
end