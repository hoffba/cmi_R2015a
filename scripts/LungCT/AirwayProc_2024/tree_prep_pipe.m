function [B,N,L_surfs] = tree_prep_pipe(caseID,fn_seg,fn_airways,outD,voxsz)
% Inputs:
%      caseID       Case ID
%      fn_seg       Lobe segmentation
%      fn_airways   Airways map
%      outD         Case directory for saving results
%      voxsz        Voxel dimensions
% Outputs:
%      B            Branches, and branch-related stats [Node1, Node2, Radius, Horsfield, Strahler, Generation]
%      N            Nodes [ID, Xi, Yi, Zi]
%      lobe_surf    Lobe surface triangulations

warning('off','MATLAB:triangulation:PtsNotInTriWarnId')
viewdir = [1,0,0];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 0. load data (ariways and lobes)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Lobe Segmentation
if ischar(fn_seg) && isfile(fn_seg)
    L = double(niftiread(fn_seg));
else
    L = fn_seg;
end
dim = size(L);

% Airways
if ischar(fn_airways) && isfile(fn_airways)
    [A,~,fov,~,~] = logical(readNIFTI(fn_airways));
    voxsz = fov./size(A,1:3);
elseif nargin==5
    A = logical(fn_airways);
else
    error('Invalid inputs\n');
end

% Check that marix size matches
if any(size(A)~=dim)
    error('Matrix dimensions do not match.');
end

% Check lobe IDs
lobe = getLobeTags(L);
lobe = lobe(cellfun(@(x)isscalar(x),{lobe.val})); % Use singular regions

% Airway mask
A = logical(A);
[A_D,A_D_idx] = bwdist(~A,'euclidean'); % to use for radii approximation

% File names for saved data:
fn_realtree = fullfile(outD,'RealTree.mat');
if isfile(fn_realtree)
    p = load(fn_realtree);
    B = p.B;
    N = p.N;
    clear p;
    nN = size(N,1);
    % nB = size(B,1);
else

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % ~~ major airway centreline graph
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    A_vox = im2coords3Dbin(A);
    A_cl = bwskel(A);
    A_cl_vox = get_M(A_cl,A_D,A_D);
    A_cl_vox = A_cl_vox(:,1:4);
    
    hf = figure('Name','Centreline Voxels'); ha = axes(hf);
    plot3(ha,A_cl_vox(:,1),A_cl_vox(:,2),A_cl_vox(:,3),'.');
    axis(ha,'equal');
    view(ha,viewdir);
    saveas(hf,fullfile(outD,[caseID,'.cl_vox.fig']));
    close(hf)
    
    [~,node,link] = Skel2Graph3D(A_cl,0);
    
    hf = figure('Name','Centreline Graph'); ha = axes(hf);
    plot_graph(ha,node,link)
    A_bound = boundary(A_vox,1); % N.B. not really tight enough, but low priority
    A_tri = triangulation(A_bound,A_vox);
    trisurf(A_tri,'FaceColor','black','FaceAlpha',0.1,'EdgeColor','none');
    grid(ha,'on');
    axis(ha,'equal');
    view(ha,viewdir);
    saveas(hf,fullfile(outD,[caseID,'.cl_graph.fig']));
    close(hf);
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % 2. formatting data into branches and nodes, more familiar structures
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Initialize branches output (B)
    nB = size(link,2);
    B = zeros(nB,3); % [Node1, Node2, Radius]
    B(:,1) = [link.n1]';
    B(:,2) = [link.n2]';
    
    % Initialize nodes output (N)
    nN = numel(node);
    N = [ (1:nN); node.comx; node.comy; node.comz ]';

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % 3. graph link radii assignment
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    % Calculate branch radii
    d2 = @(v1,v2) sqrt(sum((v1-v2).^2,2));
    for i = 1:nB
        brc = link(i);
        pnts = brc.point;
        len = length(pnts);
        
        c = (len+1)/2;
        lim = round([c-len/4 c+len/4]);
        pnts = pnts(lim(1):lim(2));
        
        tmp = zeros(1,length(pnts));
        
        for j = 1:length(tmp)
            [x,y,z] = ind2sub(dim,pnts(j));
            CLpnt = voxsz.*[x,y,z]; % point on centreline
            ind = A_D_idx(x,y,z);
            [x,y,z] = ind2sub(dim,ind);
            NRpnt = voxsz.*[x,y,z]; % nearest point from bwdist
            tmp(j) = d2(CLpnt,NRpnt);                   % Euclidean distance in mm between points
        end
        
        B(i,3) = mean(tmp); % approximate radius in mm from mean of sampled midpoints
    end
    
    % Save centerline data
    B_label = {'ID','N_Prox','N_Dist','Radius'};
    save(fn_realtree,'N','B','B_label')

end

%% 4. tree orientation

% 1. distal nodes
distN = N(histcounts(B(:,1:2),nN)==1,1);

% 2. distal branches
distB = sum(ismember(B(:,1:2),distN),2)==1;
idT = find(B(:,3)==max(B(distB,3)) & distB,1); % root branch (assuming trachea)
B = [B(idT,:) ; B(1:idT-1,:) ; B(idT+1:end,:)]; % place root first in B

% 3. check root is oriented proximal to distal
if nnz(B(:,1:2)==B(1,2))==1 % check if external node is 'distal' here
    B(1,1:2) = [B(1,2) B(1,1)]; % transpose if required
end

% 4. orient tree relative to root branch
G = graph(B(:,1),B(:,2));
TR = shortestpathtree(G,B(1,1)); %e=TR.Edges;
B2 = TR.Edges.EndNodes;

hf = figure; ha = axes(hf); view(viewdir); hold(ha,"on");
for i = 1:size(B2,1)
    try
        pnt = [N(B2(i,1),2) N(B2(i,1),3) N(B2(i,1),4)];
        vec = [N(B2(i,2),2)-N(B2(i,1),2) N(B2(i,2),3)-N(B2(i,1),3) ...
               N(B2(i,2),4)-N(B2(i,1),4)];
        quiver3(ha,pnt(1),pnt(2),pnt(3),vec(1),vec(2),vec(3),'b-',...
            'LineWidth',2,'MaxHeadSize',0.5);
    catch err
        disp(getReport(err))
        continue
    end
end
saveas(hf,fullfile(outD,[caseID,'.oriented_tree.fig']));
close(hf);

B2 = flip(B2,1); % place root at top, tree is now oriented

% 5. include radii from B in oriented tree
% Flip nodes if needed
[a,b] = ismember(B2,B(:,1:2),'rows');
B2(a,3) = B(b(logical(b)),3);
[~,b] = ismember(B2(:,[2,1]),B(:,1:2),'rows');
B2(~a,3) = B(b(logical(b)),3);
B = B2;
nB = size(B,1);
% complete conversion, B is now an oriented branch matrix with NOD as
% associated matrix of nodes, with the root (trachea) as first element

%% obtaining lobe surface maps in mm (to match tree)
%
nl = numel(lobe);
L_surfs = cell(1,nl);
for i = 1:nl

    [xx,yy,zz] = ind2sub(dim,find(L==lobe(i).val));
    tmp = [xx,yy,zz].*voxsz; clear xx yy zz
    
    s = size(tmp);
    tmp2 = datasample(tmp,floor(s(1)/10)); % random sample of 10%
    tmp3 = boundary(tmp2); % obtain indices for boundary set of voxels
    TR = triangulation(tmp3,double(tmp2));
    L_surfs{i} = TR;
end
clear tmp


%% tree generation, Strahler and Horsfield order assignment
%

H = calc_ho(B(:,1:2)); B(:,6) = H(:,end);
S = calc_so(B(:,1:2)); B(:,5) = S(:,end);
% calculating generation no. given a recursive function defined in footer
for i = 1:nB
    B(i,4) = gen_no(B(i,1:2),B(:,1:2));
end

B_label = {'Node1','Node2','Radius','Generation','Strahler','Horsfield'};
hf = figure; t = tiledlayout(hf,1,3);
for i = 4:6
    plot_tree(nexttile(t),B,B_label{i},N,i,'Raw Tree Orders',L_surfs);
end
saveas(hf,fullfile(outD,[caseID,'.raw_tree_orders.fig']));
close(hf);

hf = figure; t = tiledlayout(hf,1,3);
T = {'Generation' 'Strahler' 'Horsfield'};
for i = 4:6
    ha = nexttile(t);
    h = histogram(ha,B(:,i));
    h.BinWidth = 1;
    title(ha,T{i-3})
end
saveas(hf,fullfile(outD,[caseID,'.raw_tree_hists.fig']));
close(hf);

%% convert nodes to mm and calculate branch lengths
%
% 0. add branch lengths to end (7) of B
%
% first converting NOD to mm
N(:,2:4) = N(:,2:4).*voxsz;

% now branch length assignment using d2 metric
B(:,7) = d2( N(B(:,1),2:4) , N(B(:,2),2:4) );

hf = figure; t = tiledlayout(hf,1,2);
ha = nexttile(t);
histogram(ha,B(:,3));
title(ha,'Radii')
xlabel(ha,'airway radius (mm)')
ha = nexttile(t);
histogram(ha,B(:,7));
title(ha,'Lengths')
xlabel(ha,'airway length (mm)')
saveas(hf,fullfile(outD,[caseID,'.pre_QC_lengths_and_radii.fig']));
close(hf);

%% branch QC algorithm
%

% Idea (Brody): using thresholds t1 and t2 (set below), remove branches
% through comparison of their length (L) to parent branch length (Lp)

% L < t1*Lp => remove
% t1*Lp <= L <= t2*Lp => keep (+ match length with sibling <test>)
% L > t2*Lp => accept

% After each pass, a clean up step will follow removing interm. nodes
% Process will conintue until no terminal branches flagged for alteration

B2=B; t1=0.1; t2=0.4; flg=1;
L_fix = 2; % addition to algorithm, a fixed minimum length (mm) for all branches
% Initialize QC tracking
qc_log = [];

while flg % whole process loop
   
    n = size(B2,1); elim = true(n,1); flg = 0;
    
    % step 1: testing terminal branches to alter based on length
    for i = 2:n % always skip root
        b = B2(i,:);
        if not(ismember(b(2),B2(:,1))) % test if terminal
            B_P = B2(B2(:,2)==b(1),:); % obtain parent branch
            L_P = B_P(7); % parent branch length
            current_gen = gen_no(b(1:2), B2(:,1:2));
            
            if b(7)<t1*L_P || b(7)<L_fix
                % Log standard QC removal
                qc_log = [qc_log; current_gen, b(3), b(7), 0]; % reason: 0 = standard QC
                flg = 1; % put flag back up to continue trimming
                elim(i) = 0; % mark for elimination
            else
                % Check for terminal branches in early generations (should not exist anatomically)
                if current_gen >= 2 && current_gen <= 4
                    % Log generation 2-4 removal
                    qc_log = [qc_log; current_gen, b(3), b(7), 1]; % reason: 1 = generation 2-4
                    % Any terminal branch in generations 2-4 is suspicious - remove it
                    flg = 1;
                    elim(i) = 0; % mark for elimination
                end
            end
            if b(7) <= t2*L_P
                % take action on assumption of 'truncated branch'
            end
        end
    end
    
    B2 = B2(elim,:); % remove marked branches
    
    % cleaning of intermediate nodes...
    % idea check all branches except the root for number of children, if
    % that number is 1, merge parent and child into one new branch,
    % obtaining new properties as averages of radii and addition of lengths
    
    flg_clean = 1;
    while flg_clean 
        n = size(B2,1); flg_clean = 0;
        for i = 2:n
            j = B2(:,1)==B2(i,2);
            if nnz(j)==1 % branch has one child
                % i and ind now index rows of B2 to replace with merged brch
                B2(i,2) = B2(j,2);
                B2(i,3) = (B2(i,3)+B2(j,3))/2;
                B2(i,7) = B2(i,7)+B2(j,7);
                B2(j,:) = [];
                flg_clean = 1;
                break
            end
        end
    end
end
B = B2;
% re-evaluate Strahler and Horsfield orders
H = calc_ho(B(:,1:2)); B(:,6) = H(:,end);
S = calc_so(B(:,1:2)); B(:,5) = S(:,end);
nB = size(B,1);
for i = 1:nB
    B(i,4) = gen_no(B(i,1:2),B(:,1:2));
end
% remove redundant nodes from NOD?
% Leaving in for now to maintain row ID agreement


%% final visual check
%
hf = figure('Name','Final Visual'); ha = axes(hf);
plot_tree(ha,B,B_label,N,4,'Final Visual',L_surfs);
for i = 1:nl
    trisurf(L_surfs{i},'FaceColor','b','FaceAlpha',0.05,'EdgeColor','none');
end
saveas(hf,fullfile(outD,[caseID,'.final_visual.fig']));
close(hf);


%% saving data to pass to Brody / use for airway growing algorithm
%
% basic approximation given by multiplying by average of pixel dims, n.b.
% typically z dim is different to x and y dims which themselves match

lobe_ids = [lobe.val];
save(fullfile(outD,[caseID,'.tree_data.mat']),'B','N','lobe_ids','L_surfs','voxsz');

% Save QC report
if ~isempty(qc_log)
    qc_table = array2table(qc_log, 'VariableNames', {'Generation', 'Radius_mm', 'Length_mm', 'Removal_Type'});
    qc_table.Standard_QC = qc_table.Removal_Type == 0;
    qc_table.Gen_2to4_Removal = qc_table.Removal_Type == 1;
    writetable(qc_table, fullfile(outD, [caseID, '_QC_report.csv']));
    fprintf('QC Report: %d standard removals, %d gen 2-4 removals\n', ...
        sum(qc_table.Standard_QC), sum(qc_table.Gen_2to4_Removal));
end

warning('on','MATLAB:triangulation:PtsNotInTriWarnId')

% generation no., given all roads lead to root (B(1,:)).
function y = gen_no(b,B) % given branch b from B, output generation no.
    if b(1,1) == B(1,1) % root, define as generation 1
        y = 1;
    else % recursion using parent branch
        b = B( B(:,2)==b(1,1) ,:);
        y = 1 + gen_no(b,B);
    end