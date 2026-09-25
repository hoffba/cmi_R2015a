function iowa_preprocess2(datapath,procdir)
% Takes airways data from Iowa and saves it in the format required by MiTAP

try

    [~,ID] = fileparts(datapath);

    fprintf('Preprocessing %s ... \n',ID)

    if ~isfolder(procdir)
        mkdir(procdir);
    end

    airdir = fullfile(procdir,[ID,'_.AirwayProc']);
    if ~isfolder(airdir)
        mkdir(airdir);
    end

    % Check that all necessary files are there
    fn_B =          fullfile(datapath,[ID,'_Branches.csv']);
    fn_N =          fullfile(datapath,[ID,'_Junctions.csv']);
    fn_airways =    fullfile(datapath,'Airway_Mask.nii.gz');
    fn_CL =         fullfile(datapath,'Airway_Centerline.nii.gz');

    % File to save branches and nodes to
    fn_realtree = fullfile(airdir,'RealTree.mat');

    if isfile(fn_B) && isfile(fn_N) && isfile(fn_CL) && isfile(fn_airways)

        fprintf('    Compiling branches and nodes ...\n');

        [A,~,fov,~,~] = readNIFTI(fn_airways);
        d = size(A);
        voxsz = fov./d;

        % Branches:
        B = table2array(readtable(fn_B));
        B(:,1) = []; % Don't need a branch ID
        nB = size(B,1);

        % Nodes - need to convert to MiTAP xyz coordinates
        % ** Assumes all nodes are on centerline **
        N = readtable(fn_N);
        N = [N.Node,N.X+1,N.Y+1,N.Z+1];

        % Prep skeleton for radius calculations
        CL = logical(readNIFTI(fn_CL));
        [~, skel_D_idx] = bwdist(CL);

        % Generate graph of airway tree and calculate radius
        % Initialize skeleton structure
        % skel = struct('pStart',  mat2cell(N(B(:,1),2:4),ones(1,nB),3),...% Node location at start of branch
        %               'pEnd',    mat2cell(N(B(:,2),2:4),ones(1,nB),3),...% Node location at end of branch
        %               'points',  repmat({[]},nB,1),...                   % Centerline points of branch
        %               'gen',     repmat({nan},nB,1),...                  % Generation of branch
        %               'parent',  repmat({[]},nB,1),...                   % Parent branch
        %               'terminal',repmat({false},nB,1),...                % Flag indicating a terminal branch
        %               'complete',repmat({false},nB,1),...                % Flag for completed branches
        %               'active',  repmat({false},nB,1));                  % Flag indicating branches ready to process 
        B_term = false(nB,1);   % Flag for terminal branch
        B_gen = zeros(nB,1);    % Branch generation
        B_points = cell(nB,1);  % Branch points

        % Start at top of trachea (highest z-coordinate)
        N_start = N(N(:,4)==max(N(:,4)),1);
        iB = find(any(B==N_start,2),1); % Current branch index
        B_gen(iB) = 1;
        CLtemp = CL; % We will zero the CL as we go

        % Remove starting point from CL
        p = N(N_start,2:4);
        CLtemp(p(1),p(2),p(3)) = 0;

        hf = figure;
        im = imshow(flipud(squeeze(any(CL,1))')); ha = im.Parent;

        Bgrow = struct('StartNode',N_start,'parent',nan,'points',p);
        while ~isempty(Bgrow)

            imshow(flipud(squeeze(any(CL&~CLtemp,1))'),'Parent',ha);
            pause(0.01)

            % Growing from first branch until empty
            p = Bgrow(1).points(end,:);

            % Check if end of branch
            N_end = N(ismember(N(:,2:4),p,'rows'),1);

            % Look in window around current point
            r = 2;
            w = CLtemp(max(p(1)-r,1):min(p(1)+r,d(1)),...
                       max(p(2)-r,1):min(p(2)+r,d(2)),...
                       max(p(3)-r,1):min(p(3)+r,d(3)));
            ind = find(w);
            ngrow = numel(ind);

            if ngrow
                % Translate to image indices
                [ii,jj,kk] = ind2sub(size(w),ind);
                p_next = p+[ii,jj,kk] - (r+1);

                if ngrow>1 && isempty(N_end) % Tree branches, but this is not yet the end node
                    % Check next points for end node
                    [~,Locb] = ismember(N(:,2:4),p_next,'rows');
                    if any(Locb)
                        % Only grow to this end node
                        p_next = p_next(Locb(Locb>0),:);
                        ngrow = 1;
                    else
                        disp('Error')
                    end
                end

                % Remove points from CL
                CLtemp(p_next(:,1),p_next(:,2),p_next(:,3)) = 0;
            end
            if ngrow==1
                Bgrow(1).points(end+1,:) = p_next;
            else % End of current branch

                % Find node index
                N_end = find(ismember(N(:,2:4),p,'rows'),1);
                if isempty(N_end)
                    disp('Error')
                else
                    N_end = N(N_end,1);
                end
                % Set final points for current branch
                B_ind = find(ismember(B,[Bgrow(1).StartNode,N_end]),1);
                B_points{B_ind} = Bgrow(1).points(2:end,:);
                if ~isnan(Bgrow(1).parent)
                    B_gen(B_ind) = B_gen(Bgrow(1).parent) + 1;
                end

                if ngrow>1
                    % Add new branches
                    for i = 1:ngrow
                        Bgrow(end+1) = struct('StartNode',N_end,'parent',B_ind,'points',p_next(i,:)); %#ok<AGROW>
                    end
                else % Terminal
                    B_term(B_ind) = true;
                end

                % Remove current branch from growing list
                Bgrow(1) = [];
            end
        end

    end
catch err
    disp('break')
end








%         % Quick graph of skeleton
%         [~,node,link] = Skel2Graph3D(CL,0);
%         link_point = {link.point}';
%         link_point_all_idx = [link_point{:}]';
%         [rr,cc,ss] = ind2sub(d,link_point_all_idx);
%         link_point_all = [rr,cc,ss];
% 
%         % Make sure all nodes are on skeleton graph
%         missing_point = find(~ismember(N_idx,link_point_all_idx));
%         k = dsearchn(link_point_all,N(missing_point,2:4));
%         for ik = 1:length(k)
% 
%             % Find link associated with nearest point
%             k_idx = link_point_all_idx(k(ik));
%             klink = find(cellfun(@(x)ismember(k_idx,x),link_point),1);
% 
%             % Find location relative to link points list
%             point_list_idx = link_point{klink}';
%             np = length(point_list_idx);
%             [rr,cc,ss] = ind2sub(d,point_list_idx);
%             point_list = [rr,cc,ss];
%             [~,dist] = dsearchn(N(missing_point(ik),2:4),point_list);
%             pdiff = sqrt(sum(diff(point_list,1).^2,2));
%             kpt = find(dist==min(dist));
%             mp_idx = N_idx(missing_point(ik));
%             % Decide where to put it
%             if kpt == 1
%                 if dist(2)<pdiff(1)
%                     btwn = [1 2];
%                 else
%                     btwn = [0 1];
%                 end
%             elseif kpt == numel(dist)
%                 if dist(end-1)<pdiff(end)
%                     btwn = [np-1 np];
%                 else
%                     btwn = [np np+1];
%                 end
%             else
%                 kmin = dist==min(dist);
%                 knext = find(dist==min(dist(~kmin)),1);
%                 kmin = find(kmin,1);
%                 btwn = [min(kmin,knext),max(kmin,knext)];
%             end
%             link_point{klink} = [ point_list_idx(1:btwn(1)) ; mp_idx ; point_list_idx(btwn(2):end) ]';
% 
%         end
% 
%         % Find voxels just outside airways for radius calculation
%         A_border_idx = find(imdilate(A,strel('sphere',1)) & ~A);
%         [rr,cc,ss] = ind2sub(d,A_border_idx);
%         A_border_point = [rr,cc,ss];
% 
%         % Loop over branches and calculate radii
%         fprintf('    Calculating branch mean radii ...\n')
%         for iB = 1:nB
% 
%             % start node
%             ind = N(:,1)==B(iB,1);
%             start_node = N(ind,2:4);
%             start_node_idx = N_idx(ind);
%             % Check if node is on the skeleton
%             if ~CL(start_node(1), start_node(2), start_node(3))
%                 idx = sub2ind(d, start_node(1), start_node(2), start_node(3));
%                 start_node_idx = skel_D_idx(idx);
%                 % [start_node(1), start_node(2), start_node(3)] = ind2sub(d, start_node_idx);
%             end
% 
%             % end node
%             ind = N(:,1)==B(iB,2);
%             end_node =   N(ind,2:4);
%             end_node_idx = N_idx(ind);
%             % Check if node is on the skeleton
%             if ~CL(end_node(1), end_node(2), end_node(3))
%                 idx = sub2ind(d, end_node(1), end_node(2), end_node(3));
%                 end_node_idx = skel_D_idx(idx);
%             end
% 
%             % Find SkelGraph links containing start and end nodes
%             ind = cellfun(@(x)any(ismember([start_node_idx,end_node_idx],x)),link_point);
% 
%             % Create graph with link points
%             gpoint = unique([link_point{ind}]);
%             ngp = length(gpoint);
%             G = graph;
%             G = addnode(G,ngp);
%             for igp = 1:ngp
%                 [r, c, p] = ind2sub(d, gpoint(igp));
% 
%                 % Use 26-connectivity for neighbors.
%                 [neighbor_r, neighbor_c, neighbor_p] = ndgrid(r-1:r+1, c-1:c+1, p-1:p+1);
%                 neighbors_idx = sub2ind(d, neighbor_r(:), neighbor_c(:), neighbor_p(:));
%                 neighbors_idx(14) = [];  % Exclude point itself
% 
%                 % Find neighbors in set of skeleton points
%                 ind = find(ismember(gpoint,neighbors_idx));
% 
%                 for i_neighbor = ind'
%                     % Add edge only once to avoid duplicates.
%                     G = addedge(G, igp, i_neighbor);
%                 end
%             end
% 
%             % Find the shortest path between start and end nodes
%             % Find the graph node IDs for the start and end points.
%             start_id = find(start_node_idx==gpoint);
%             end_id = find(end_node_idx==gpoint);
%             path_nodes = shortestpath(G, start_id, end_id);
% 
%             % Convert path back to voxel coordinates
%             path_linear_indices = gpoint(path_nodes);
%             [path_rows, path_cols, path_planes] = ind2sub(d, path_linear_indices');
%             path_points = [path_rows, path_cols, path_planes];
% 
%             % Calculate mean radius of branch
%             % - only calculate middle section of branch
%             if ~isempty(path_points)
% 
%                 % Use points with margin away from nodes
%                 np = size(path_points,1);
%                 c = (np+1)/2;
%                 lim = round([c-np/4 c+np/4]);
%                 path_points = path_points(lim(1):lim(2),:);
%                 % np = size(path_points,1);
% 
%                 % Find distance to nearest airway border
%                 [~,dist] = dsearchn(A_border_point,path_points);
% 
%                 B(iB,3) = mean(dist); % approximate radius in mm from mean of sampled midpoints
%             else
%                 B(iB,3) = nan;
%             end
%         end
% 
%         % Save tree to file for MiTAP to load
%         B_label = {'N_Prox','N_Dist','Radius'};
%         save(fn_realtree,'N','B','B_label');
%     end
% 
%     % Copy certain files to the processing directory
%     fn_copy = {'CT_Expiratory',         'exp';...
%                'CT_Inspiratory',        'ins';...
%                'Lung_Lobes_Expiratory', 'exp.label';...
%                'Lung_Lobes_Inspiratory','ins.label';...
%                'Airway_Mask',           'airways'};
%     for j = 1:size(fn_copy,1)
%         fprintf('    Copying file: %s to %s\n',fn_copy{j,:});
%         fn = fullfile(datapath,ID{i},[fn_copy{j,1},'.nii.gz']);
%         if isfile(fn)
%             copyfile(fn,fullfile(procdir,[ID{i},'_.',fn_copy{j,2},'.nii.gz']),CopyLinkBehavior="resolve");
%         end
%     end
% 
%     fprintf('  done\n')
% 
% catch err
%     disp('breakpoint')
% end
% 
% 
% 
% 
% 
