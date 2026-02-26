function iowa_preprocess(datapath)
% Takes airways data from Iowa and saves it in the format required by MiTAP


[~,ID] = fileparts(datapath);

fprintf('Preprocessing %s ... \n',ID)

airdir = fullfile(datapath,[ID,'.AirwayProc']);
if ~isfolder(airdir)
    mkdir(airdir);
end

% Check that all necessary files are there
fn_B =          fullfile(datapath,[ID,'_Branches.csv']);
fn_N =          fullfile(datapath,[ID,'_Junctions.csv']);
fn_airways =    fullfile(datapath,'Airway_Mask.nii.gz');
fn_CL =         fullfile(datapath,'Airway_Centerline.nii.gz');
fn_seg =        fullfile(datapath,'Lung_Lobes_Inspiratory.nii.gz');

if isfile(fn_B) && isfile(fn_N) && isfile(fn_CL) && isfile(fn_airways) && isfile(fn_seg)

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
    N_ind = sub2ind(d,N(:,2),N(:,3),N(:,4));

    % Re-save skeleton for pipeline processing
    copyfile(fn_CL,fullfile(airdir,[ID,'.Centerline.nii.gz']))
    CL = logical(readNIFTI(fn_CL));
    
    % Find CL points belonging to each branch    
    G = bwgraph(CL);
    points = cell(nB,1);
    fprintf('Processing branches (%d) ... 000',nB);
    for i = 1:nB
        fprintf('\b\b\b%03d',i)
        % Find branch CL points
        path_points = shortestpath(G,N_ind(B(i,1)),N_ind(B(i,2)));
        points{i} = path_points(2:end-1); % Exclude node points
    end
    fprintf('\n');

    % Find missing branches
    %  ** BH found that one branch was not included in the branch list for all cases
    missingNi = find(~ismember(N(:,1),B)); % Nodes not in the Branches list
    nmiss = numel(missingNi);
    B(end+(1:nmiss),:) = nan;
    points(end+(1:nmiss)) = {[]};
    for i = 1:nmiss
        [TR,D] = shortestpathtree(G,N_ind(missingNi(i)),N_ind,'OutputForm','cell');
        D(D==0) = nan;
        [~,nearestTargetIdx] = min(D,[],'omitnan');
        path_points = TR{nearestTargetIdx};
        B(nB+i,:) = [nearestTargetIdx,N(missingNi,1)];
        points{nB+i} = path_points;
    end

    % Put tracheal branch first
    nB = size(B,1);
    activeN = N(find(N(:,4)==max(N(:,4)),1),1);
    Bi = find(any(B(:,1:2)==activeN,2),1);
    ind = [Bi,1:Bi-1,Bi+1:nB];
    B = B(ind,:);
    points = points(ind);

    % Save tree data
    N(:,2:4) = N(:,2:4).*voxsz;
    save(fullfile(airdir,[ID,'.RealTree.mat']),'B','N','points','d','voxsz')

    % Save QC figure
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
    plot3(ha,CL_bk(:,1),CL_bk(:,2),CL_bk(:,3),'k.');
    plot3(ha,CL_rd(:,1),CL_rd(:,2),CL_rd(:,3),'r.');
    saveas(hf,fullfile(datapath,[ID,'_AirwayTree.tif']));
    close(hf);

end

