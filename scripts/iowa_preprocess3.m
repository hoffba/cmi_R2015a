function iowa_preprocess3(datapath,procdir)
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
        A = logical(A);
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
        % [~, skel_D_idx] = bwdist(CL);

        % Find voxels just outside airways for radius calculation
        A_border_point = find(imdilate(A,strel('sphere',1)) & ~A);
        [rr,cc,ss] = ind2sub(d,A_border_point);
        A_border_point = [rr,cc,ss] .* voxsz;

        % Generate graph of airway tree and calculate radius
        G = bwgraph(CL);

        % Find CL points belonging to each branch
        points = cell(nB,1);
        for i = 1:nB

            fprintf('%03d\n',i)

            % Find branch CL points
            BN = N(B(i,1:2),2:4);
            BN = sub2ind(d,BN(:,1),BN(:,2),BN(:,3));
            path_points = shortestpath(G,BN(1,:),BN(2,:));

            % Calculate average radius
            %   Use points with margin away from end nodes
            np = numel(path_points);
            c = (np+1)/2;
            lim = round([c-np/4 c+np/4]);
            mean_points = path_points(lim(1):lim(2));
            [rr,cc,ss] = ind2sub(d,mean_points');

            % Find distance to nearest airway border
            [~,dist] = dsearchn(A_border_point,[rr,cc,ss] .* voxsz);

            B(i,3) = mean(dist); % approximate radius in mm from mean of sampled midpoints

            points{i} = path_points;
        end

    end

    CL_tree = false(d);
    CL_tree([points{:}]) = true;
    [rr,cc,ss] = ind2sub(d,[points{:}]');
    CL_bk = [rr,cc,ss] .* voxsz;
    [rr,cc,ss] = ind2sub(d,find(CL & ~CL_tree));
    CL_rd = [rr,cc,ss] .* voxsz;

    % Prep figure / axes
    hf = figure('Name','Airway Tree'); ha = axes(hf);
    axis(ha,'equal');
    hold(ha,'on');
    grid(ha,'on');
    view(ha,[1,0,0]);
    title(ha,'Iowa Proc');
    plot3(ha,CL_bk(:,1),CL_bk(:,2),CL_bk(:,3),'k.')
    plot3(ha,CL_rd(:,1),CL_rd(:,2),CL_rd(:,3),'r.')

    disp('Check')
catch err
    disp('Error')
end