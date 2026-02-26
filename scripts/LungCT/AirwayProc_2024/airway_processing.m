% airway_processing.m
function [B,N,L_surfs] = airway_processing(pID,seg,A,voxsz,procdir)
%% function airway_processing(pID,ins,seg,A,voxsz,img,fcn,procdir)
%   Inputs: pID =       string ID of subject/case (base filename)
%           seg =       inspiratory lung segmentation
%           A =         airways segmentation
%           voxsz =     voxel dimensions
%           procdir =   directory to save results to
%   Outputs: B = [ID, N_prox, N_dist, Radius, Gen, Strahl, Hors, Lobe, Tag]
%            N = [ID, x, y, z]
%            lobe_surf = structure containing lobe surface triangulations
%
%   Tag values:
%       0 = Real non-terminal branch
%      <0 = Real terminal branch (negative index = terminal ID)
%      >0 = Simulated branch (positive value = which terminal it belongs to)

viewdir = [1,0,0];

if ischar(seg)
    seg = cmi_load(1,[],seg);
end
if ischar(A)
    A = cmi_load(1,[],A);
end

% Check that inputs are the same size
dim = size(seg,1:3);
if ~all(size(A,1:3) == dim)
    error('All input dimensions must match.')
end
fov = dim.*voxsz;

%% Create folder for outputs
outD = fullfile(procdir,[pID,'.AirwayProc']);
if ~isfolder(outD)
    mkdir(outD);
end

%% Define file names
fn = struct('real',       fullfile(outD,[pID,'.RealTree.mat'])     ,...
            'sim',        fullfile(outD,[pID,'.SimTree.mat'])      ,...
            'lobe',       fullfile(outD,[pID,'.LobeSurf.mat'])     ,...
            'centerline', fullfile(outD,[pID,'.Centerline.nii.gz']));
qc = struct('cl_vox',     fullfile(outD,[pID,'_cl_vox.tif'])        ,...
            'cl_graph',   fullfile(outD,[pID,'_cl_graph.tif'])      );

%% Process lobes (for QC figure surface rendering)
if isfile(fn.lobe)
    fprintf('Reading lobe surfaces from file\n');
    p = load(fn.lobe);
    L_surfs = p.L_surfs;
    % lobe_ids = p.lobe_ids;
else
    fprintf('Generating lobe surfaces ... ');
    t = tic;
    [L_surfs,lobe_ids] = getLobeSurfs(seg,voxsz);
    save(fn.lobe,'lobe_ids','L_surfs','voxsz');
    t = toc(t);
    fprintf('done (%d:%02d)\n',floor(t/60),rem(t,60));
end

%% Skeletonize the airways
if isfile(fn.centerline)
    fprintf('Reading centerline from file: %s\n',fn.centerline);
    CL = readNIFTI(fn.centerline);
    CL = logical(CL);
else
    fprintf('Finding centerline from airways segmentation ... ');
    t = tic;
    CL = findCenterline(logical(A));
    saveNIFTI(fn.centerline,CL,{'Centerline'},fov,eye(4))
    t = toc(t);
    fprintf('done (%d:%02d)\n',floor(t/60),rem(t,60));
end

%% QC figure shows centerline
[xx,yy,zz] = ind2sub(dim,find(CL));
A_cl_vox = [xx,yy,zz] .* voxsz;
hf = figure('Name','Centreline Voxels'); ha = axes(hf);
plot3(ha,A_cl_vox(:,1),A_cl_vox(:,2),A_cl_vox(:,3),'.');
axis(ha,'equal');
view(ha,viewdir);
saveas(hf,qc.cl_vox);
close(hf)

%% Convert skeleton to nodes and branches
if isfile(fn.real)
    fprintf('Reading airway tree from file: %s\n',fn.real)
    p = load(fn.real);
    B = p.B;            % Branches described by Node IDs
    N = p.N;            % Node positions in image-xyz (mm)
    points = p.points;  % Branch centerline point indices (not including end nodes)
else
    fprintf('Generating airway tree ... ');
    t = tic;
    [B,N,points] = skel2tree(CL,voxsz);
    save(fn_real,'B','N','points');
    t = toc(t);
    fprintf('done (%d:%02d)\n',floor(t/60),rem(t,60));
end

%% QC figure shows tree in airway rendering
% Generate airway surface
A_surf = getLobeSurfs(A,voxsz);
% Generate figure
hf = figure('Name','Centreline Graph'); ha = axes(hf);
plot_tree(ha,B,N,[],A_surf);
saveas(hf,fullfile(outD,[pID,'.cl_graph.fig']));
close(hf);

%% QC rearrange and fix tree errors

% Remove branches with 0 length
ind = B(:,1)==B(:,2);
if nnz(ind)
    fprintf('Removing %d branches with zero length.\n',nnz(ind));
    B(ind,:) = [];
end

% Place tracheal branch first in the list
%  * slices go from 0 (bottom) to N (top) of lungs
Nroot = N(find(N(:,4)==max(N(:,4)),1),1);
Bi = find(any(B(:,1:2)==Nroot,2),1);
B = [B(Bi,:);B(1:Bi-1,:);B(Bi+1:end,:)];
if B(1,1)~=Nroot
    B(1,:) = B(1,[2,1]);
end

% Orient tree branches to be proximal-distal
fprintf('Orienting branches starting at trachea ... ');
t = tic;
G = graph(B(:,1),B(:,2));
TR = shortestpathtree(G,B(1,1));
B = TR.Edges.EndNodes;
t = toc(t);
fprintf('done (%d:%02d)\n',floor(t/60),rem(t,60));

% Combine branches that don't split
[uval, ~, ind] = unique(B(:));          % 1. Find the unique values and their positions
counts = histcounts(ind, 1:max(ind)+1); % 2. Count the occurrences of each unique value
midnode = uval(counts == 2);            % 3. Find values that appeared exactly twice
ncombine = numel(midnode);
if ncombine
    fprintf('Combining limbs that do not branch (%d) ...',ncombine);
    t = tic;
    for i = 1:ncombine
        B_prox = find(B(:,2)==midnode(i),1);
        B_dist = find(B(:,1)==midnode(i),1);
        B(B_prox,2) = B(B_dist,2);
        Ni = round(N(midnode,2:4)./voxsz);
        points{B_prox} = [points{B_prox} , sub2ind(dim,Ni(1),Ni(2),Ni(3)) , points{B_dist}];
        B(B_dist,:) = [];
    end
    t = toc(t);
    fprintf('done (%d:%02d)\n',floor(t/60),rem(t,60));
end

nB = size(B,1);

%% QC figure shows quiver of branch directions
hf = figure; ha = axes(hf); view(viewdir); hold(ha,"on");
for i = 1:nB
    try
        pnt = N(B(i,1),2:4);
        vec = N(B(i,2),2:4) - N(B(i,1),2:4);
        quiver3(ha,pnt(1),pnt(2),pnt(3),vec(1),vec(2),vec(3),'b-','LineWidth',2,'MaxHeadSize',0.5);
    catch err
        disp(getReport(err))
        continue
    end
end
saveas(hf,fullfile(outD,[pID,'.oriented_tree.fig']));
close(hf);


%% Calculate branch values
% Mean Radius

% Find voxels just outside airways for radius calculation
fprintf('Calculating mean branch radii ...\n')
A_border_point = find(imdilate(A,strel('sphere',1)) & ~A);
[rr,cc,ss] = ind2sub(dim,A_border_point);
A_border_point = [rr,cc,ss] .* voxsz;
R = nan(nB,1);
for i = 1:nB
    path_points = points{i};
    np = numel(path_points);
    if np % Can't calculate if there are no points
        % Use points with margin away from end nodes
        c = (np+1)/2;
        lim = round([c-np/4 c+np/4]);
        mean_points = path_points(lim(1):lim(2));
        [rr,cc,ss] = ind2sub(dim,mean_points');
    
        % Find distance to nearest airway border
        [~,dist] = dsearchn(A_border_point,[rr,cc,ss] .* voxsz);
        R(i) = mean(dist); % approximate radius in mm from mean of sampled midpoints
    end
end

% Horsfield Order
H = calc_ho(B);
% Strahler Order
S = calc_so(B);
% Generation
G = calc_gen(B);

% Concatenate into Branch matrix
B = [B,R,H,S,G];


%% Ali Namvar, 2026-02-17 - Save real terminal distal node xyz before simulation
real_term_idx = findTerminalBranches(B);
Nr_term_xyz = N(B(real_term_idx, 2), 2:4);
fprintf('Real terminal branches before simulation: %d\n', nnz(real_term_idx));

%% Simulation
t = tic;
[Bsim,Nsim,~,~] = CreateConductingZone( B(:,1:3), N, L_surfs , outD );  % Ali Namvar, 2026-02-17 - removed num_real_branches output
t = round(toc(t));
fprintf('Simulation = %d:%02d\n',floor(t/60),rem(t,60));
% Bsim ~ [ID, Node1, Node2, Radius, Gen, Strahler, Horsfield, Lobe]

%% Ali Namvar, 2026-02-17 - Tag branches using xyz coordinate matching
t = tic;
nB = size(Bsim, 1);
tag = zeros(nB, 1);

% Step 1: Find real terminal nodes in post-simulation tree by xyz matching
Ns_term = Nsim(ismember(Nsim(:,2:4), Nr_term_xyz, 'rows'), 1);

if isempty(Ns_term)
    warning('No real terminal nodes found by xyz matching - all tags will be 0');
end

% Step 2: Find branches whose distal node is a real terminal node
Br_term = find(ismember(Bsim(:,3), Ns_term));
nterm = numel(Br_term);

fprintf('\nTagging branches...\n');
fprintf('  Total branches: %d\n', nB);
fprintf('  Real terminals matched by xyz: %d / %d\n', nterm, length(real_term_idx));

% Step 3: Tag real terminals with negative numbers
tag(Br_term) = -(1:nterm);

% Step 4: Tag simulated descendants with positive numbers
for i = 1:nterm
    Bnext = ismember(Bsim(:,2), Bsim(Br_term(i),3));
    while nnz(Bnext)
        tag(Bnext) = i;
        Bnext = ismember(Bsim(:,2), Bsim(Bnext,3));
    end
end

%% Validate tagging
n_real_nonterminal = sum(tag == 0);
n_real_terminal = sum(tag < 0);
n_simulated = sum(tag > 0);

fprintf('\nTag validation:\n');
fprintf('  Tag = 0 (real non-terminal): %d\n', n_real_nonterminal);
fprintf('  Tag < 0 (real terminal): %d\n', n_real_terminal);
fprintf('  Tag > 0 (simulated): %d\n', n_simulated);
fprintf('  Unmatched terminals: %d\n', length(real_term_idx) - nterm);

if nterm < length(real_term_idx)
    warning('%d real terminals could not be matched in post-simulation tree', ...
        length(real_term_idx) - nterm);
end

t = round(toc(t));
fprintf('Branch tagging = %d:%02d\n',floor(t/60),rem(t,60));

%% Assemble output
B = [Bsim, tag];
N = Nsim;

%% Save tree
B_label = {'ID','N_Prox','N_Dist','Radius','Gen','Strahler','Horsfield','Lobe','Tag'};
dims = size(seg);

% Ali Namvar, 2026-02-17 - Updated tag_stats, removed num_real_branches
tag_stats.total_branches = nB;
tag_stats.real_nonterminal = n_real_nonterminal;
tag_stats.real_terminal = n_real_terminal;
tag_stats.simulated = n_simulated;
tag_stats.real_terminal_indices = Br_term;

save(fullfile(outD,[pID,'_AirwayTreeSim.mat']),'dims','voxsz','B','N','B_label','L_surfs','tag_stats');

fprintf('\nSaved airway tree to: %s\n', fullfile(outD,[pID,'_AirwayTreeSim.mat']));
fprintf('Final tree: %d branches (%d real, %d simulated)\n', nB, n_real_nonterminal + n_real_terminal, n_simulated);
