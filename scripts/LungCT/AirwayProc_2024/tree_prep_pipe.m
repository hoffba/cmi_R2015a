function [B,N,L_surfs] = tree_prep_pipe(caseID,fn_seg,fn_airways,outD,voxsz)
% Inputs:
%      caseID       Case ID
%      fn_seg       Lobe segmentation
%      fn_airways   Airways map
%      outD         Case directory for saving results
%      voxsz        Voxel dimensions
% Outputs:
%      B            Branches, and branch-related stats [Node1, Node2, Radius, Generation, Strahler, Horsfield, Length]
%      N            Nodes [ID, Xi, Yi, Zi]
%      lobe_surf    Lobe surface triangulations

warning('off','MATLAB:triangulation:PtsNotInTriWarnId')
viewdir = [1,0,0];

%% 0. load data (airways and lobes)

% Lobe Segmentation
if ischar(fn_seg) && isfile(fn_seg)
    L = double(niftiread(fn_seg));
else
    L = fn_seg;
end
dim = size(L);

% Airways
if ischar(fn_airways) && isfile(fn_airways)
    [A,~,fov,~,~] = readNIFTI(fn_airways);
    voxsz = fov./size(A,1:3);
elseif nargin==5
    A = fn_airways;
else
    error('Invalid inputs\n');
end

% Check that matrix size matches
if any(size(A)~=dim)
    error('Matrix dimensions do not match.');
end

% Check lobe IDs
lobe = getLobeTags(L);
lobe = lobe(cellfun(@(x)numel(x)==1,{lobe.val})); % Use singular regions

%% ============================================================
%% FIX: Select largest connected component BEFORE skeletonization
%% ============================================================
A = logical(A);

% Find all connected components in 3D
CC = bwconncomp(A, 26);  % 26-connectivity for 3D

if CC.NumObjects == 0
    error('No airway components found in the airway mask.');
elseif CC.NumObjects > 1
    % Find the largest component
    numPixels = cellfun(@numel, CC.PixelIdxList);
    [maxPixels, largestIdx] = max(numPixels);
    
    % Log information about discarded components
    totalPixels = sum(numPixels);
    discardedPixels = totalPixels - maxPixels;
    fprintf('Connected component analysis:\n');
    fprintf('  Total components found: %d\n', CC.NumObjects);
    fprintf('  Largest component: %d voxels (%.1f%% of total)\n', maxPixels, 100*maxPixels/totalPixels);
    fprintf('  Discarded: %d voxels from %d smaller components\n', discardedPixels, CC.NumObjects-1);
    
    % Create new mask with only the largest component
    A_clean = false(size(A));
    A_clean(CC.PixelIdxList{largestIdx}) = true;
    A = A_clean;
    clear A_clean;
    
    % Save QC info about component selection
    comp_info.num_components = CC.NumObjects;
    comp_info.largest_size = maxPixels;
    comp_info.total_size = totalPixels;
    comp_info.percent_kept = 100*maxPixels/totalPixels;
    comp_info.discarded_sizes = numPixels(numPixels ~= maxPixels);
    save(fullfile(outD,[caseID,'_component_QC.mat']),'comp_info');
else
    fprintf('Airway mask has single connected component - no filtering needed.\n');
end

%% 1. the major airway centreline graph

[A_D,A_D_idx] = bwdist(~A,'euclidean'); % to use for radii approximation

A_vox = im2coords3Dbin(A);
A_cl = bwskel(A, 'MinBranchLength', 10);  % Remove short spurious branches during skeletonization
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
A_bound = boundary(A_vox,1);
A_tri = triangulation(A_bound,A_vox);
trisurf(A_tri,'FaceColor','black','FaceAlpha',0.1,'EdgeColor','none');
grid(ha,'on');
axis(ha,'equal');
view(ha,viewdir);
saveas(hf,fullfile(outD,[caseID,'.cl_graph.fig']));
close(hf);


%% 2. formatting data into branches and nodes
%
% Initialize branches output (B)
nB = size(link,2);
B = zeros(nB,3); % [Node1, Node2, Radius]
B(:,1) = [link.n1]';
B(:,2) = [link.n2]';

% Initialize nodes output (N)
nN = numel(node);
N = [ (1:nN); node.comx; node.comy; node.comz ]';


%% 3. graph link radii assignment

% Calculate branch radii
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
        tmp(j) = sqrt(sum((CLpnt-NRpnt).^2,2));
    end
    
    B(i,3) = mean(tmp); % approximate radius in mm
end


%% 4. tree orientation

% 1. distal nodes
distN = N(histcounts(B(:,1:2),nN)==1,1);

% 2. distal branches
distB = sum(ismember(B(:,1:2),distN),2)==1;
idT = find(B(:,3)==max(B(distB,3)) & distB,1); % root branch (assuming trachea)
B = [B(idT,:) ; B(1:idT-1,:) ; B(idT+1:end,:)]; % place root first in B

% 3. check root is oriented proximal to distal
if nnz(B(:,1:2)==B(1,2))==1
    B(1,1:2) = [B(1,2) B(1,1)];
end

% 4. orient tree relative to root branch
G = graph(B(:,1),B(:,2));
TR = shortestpathtree(G,B(1,1));
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

B2 = flip(B2,1); % place root at top

% 5. include radii from B in oriented tree
[a,b] = ismember(B2,B(:,1:2),'rows');
B2(a,3) = B(b(logical(b)),3);
[~,b] = ismember(B2(:,[2,1]),B(:,1:2),'rows');
B2(~a,3) = B(b(logical(b)),3);
B = B2;
nB = size(B,1);

%% obtaining lobe surface maps in mm (to match tree)
nl = numel(lobe);
L_surfs = cell(1,nl);
for i = 1:nl
    [xx,yy,zz] = ind2sub(dim,find(L==lobe(i).val));
    tmp = [xx,yy,zz].*voxsz; clear xx yy zz
    
    s = size(tmp);
    tmp2 = datasample(tmp,floor(s(1)/10));
    tmp3 = boundary(tmp2);
    TR = triangulation(tmp3,double(tmp2));
    L_surfs{i} = TR;
end
clear tmp


%% tree generation, Strahler and Horsfield order assignment

H = calc_ho(B(:,1:2)); B(:,6) = H(:,end);
S = calc_so(B(:,1:2)); B(:,5) = S(:,end);
for i = 1:nB
    B(i,4) = gen_no(B(i,1:2),B(:,1:2));
end

B_label = {'Node1','Node2','Radius','Generation','Strahler','Horsfield'};
hf = figure; t = tiledlayout(hf,1,3);
for i = 4:6
    plot_tree(nexttile(t),B,N,B(:,i),L_surfs,B_label{i});
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

% First convert N to mm
N(:,2:4) = N(:,2:4).*voxsz;

% Calculate branch lengths using Euclidean distance
B(:,7) = sqrt(sum((N(B(:,1),2:4) - N(B(:,2),2:4)).^2, 2));

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

%% ============================================================
%% STANDARD BRANCH QC (Length-based)
%% ============================================================

B2 = B; 
t1 = 0.1; 
t2 = 0.4; 
flg = 1;
L_fix = 2; % minimum length (mm)
qc_log = [];

while flg
    n = size(B2,1); 
    elim = true(n,1); 
    flg = 0;
    
    for i = 2:n
        b = B2(i,:);
        if ~ismember(b(2), B2(:,1)) % is terminal
            B_P = B2(B2(:,2)==b(1),:);
            L_P = B_P(7);
            if b(7) < t1*L_P || b(7) < L_fix
                current_gen = gen_no(b(1:2), B2(:,1:2));
                qc_log = [qc_log; current_gen, b(3), b(7), 0];
                flg = 1;
                elim(i) = 0;
            end
        end
    end
    
    B2 = B2(elim,:);
    
    % Clean up single-child nodes (merge branches)
    flg_clean = 1;
    while flg_clean 
        n = size(B2,1); 
        flg_clean = 0;
        for i = 2:n
            j = B2(:,1) == B2(i,2);
            if nnz(j) == 1
                B2(i,2) = B2(j,2);
                B2(i,3) = B2(j,3); % Ali Namvar, 2026-02-17 - Use child radius (surviving airway)
                % Recalculate length from actual node positions
                B2(i,7) = sqrt(sum((N(B2(i,1),2:4) - N(B2(i,2),2:4)).^2));
                B2(j,:) = [];
                flg_clean = 1;
                break
            end
        end
    end
end
B = B2;

%% ============================================================
%% SPURIOUS EARLY TERMINAL REMOVAL
%% ============================================================
% Two-tier removal strategy:
%   Gen 1-4: Strahler_diff >= 5 ONLY (ignore length!)
%            At early generations, no real airway can be terminal
%            while sibling goes 5+ levels deeper
%   Gen 5-6: Strahler_diff >= 3 AND Length < 8mm
%            Need both criteria for smaller airways
%
% After each removal, generation numbers are recalculated!

fprintf('\n=== SPURIOUS EARLY TERMINAL QC ===\n');
fprintf('Criteria:\n');
fprintf('  Gen 1-4: Strahler_diff >= 5 (ignore length)\n');
fprintf('  Gen 5-6: Strahler_diff >= 3 AND Length < 7mm\n\n');

MAX_GEN_TO_CHECK = 6;           % Check Gen 1-6
MIN_STRAHLER_DIFF_GEN_1_4 = 5;  % For Gen 1-4: Strahler only
MIN_STRAHLER_DIFF_GEN_5_6 = 3;  % For Gen 5-6: Strahler + length
LENGTH_THRESH_GEN_5_6 = 8;      % mm - only for Gen 5-6

early_term_log = [];
flg_early = 1;
iter = 0;
max_iter = 100;  % Safety limit

while flg_early && iter < max_iter
    flg_early = 0;
    iter = iter + 1;
    n = size(B,1);
    
    % IMPORTANT: Recalculate Strahler and Generation EVERY iteration
    % This ensures correct values after any removal/merge
    S = calc_so(B(:,1:2)); 
    B(:,5) = S(:,end);
    H = calc_ho(B(:,1:2));
    B(:,6) = H(:,end);
    for i = 1:n
        B(i,4) = gen_no(B(i,1:2), B(:,1:2));
    end
    
    % Also recalculate all lengths from actual node positions
    for i = 1:n
        B(i,7) = sqrt(sum((N(B(i,1),2:4) - N(B(i,2),2:4)).^2));
    end
    
    for i = 2:n  % Skip root
        my_gen = B(i,4);
        
        % Only check early generations
        if my_gen > MAX_GEN_TO_CHECK
            continue
        end
        
        % Check if terminal (Strahler = 1)
        if B(i,5) ~= 1
            continue
        end
        
        % Verify it's actually terminal (no children)
        if ismember(B(i,2), B(:,1))
            continue
        end
        
        % Get sibling branches (branches with same proximal node)
        sibling_idx = find(B(:,1) == B(i,1) & (1:n)' ~= i);
        
        if isempty(sibling_idx)
            continue
        end
        
        % Get metrics
        sibling_strahler = max(B(sibling_idx, 5));
        my_strahler = B(i, 5);
        my_length = B(i, 7);
        
        strahler_diff = sibling_strahler - my_strahler;
        
        % Apply DIFFERENT criteria based on generation
        if my_gen <= 4
            % Gen 1-4: Strahler asymmetry ONLY (ignore length)
            % No real airway at Gen 1-4 can be terminal while sibling goes 5+ deeper
            is_spurious = strahler_diff >= MIN_STRAHLER_DIFF_GEN_1_4;
            criteria_str = sprintf('Strahler_diff >= %d', MIN_STRAHLER_DIFF_GEN_1_4);
        else
            % Gen 5-6: Need BOTH Strahler AND length criteria
            crit_strahler = strahler_diff >= MIN_STRAHLER_DIFF_GEN_5_6;
            crit_length = my_length < LENGTH_THRESH_GEN_5_6;
            is_spurious = crit_strahler && crit_length;
            criteria_str = sprintf('Strahler_diff >= %d AND Length < %.0fmm', ...
                MIN_STRAHLER_DIFF_GEN_5_6, LENGTH_THRESH_GEN_5_6);
        end
        
        % Print decision
        if is_spurious
            fprintf('REMOVING spurious terminal at Gen=%d:\n', my_gen);
        else
            fprintf('KEEPING terminal at Gen=%d (does not meet criteria):\n', my_gen);
        end
        fprintf('  Criteria for Gen %d: %s\n', my_gen, criteria_str);
        fprintf('  Strahler_diff=%d, Length=%.2f mm\n', strahler_diff, my_length);
        fprintf('  Sibling Strahler=%d, My Strahler=%d\n', sibling_strahler, my_strahler);
        if my_gen <= 4
            fprintf('  Decision: Strahler_diff(%d) >= %d? %s\n\n', ...
                strahler_diff, MIN_STRAHLER_DIFF_GEN_1_4, bool2str(is_spurious));
        else
            fprintf('  Decision: Strahler_diff(%d) >= %d? %s, Length(%.1f) < %d? %s\n\n', ...
                strahler_diff, MIN_STRAHLER_DIFF_GEN_5_6, bool2str(strahler_diff >= MIN_STRAHLER_DIFF_GEN_5_6), ...
                my_length, LENGTH_THRESH_GEN_5_6, bool2str(my_length < LENGTH_THRESH_GEN_5_6));
        end
        
        if is_spurious
            % Log: [Gen, Length, Strahler_diff, Sibling_Strahler]
            early_term_log = [early_term_log; my_gen, my_length, strahler_diff, sibling_strahler];
            
            % Remove the branch
            B(i,:) = [];
            flg_early = 1;
            break  % Restart loop - indices changed
        end
    end
    
    % Clean up single-child nodes after removal
    flg_clean = 1;
    while flg_clean 
        n = size(B,1); 
        flg_clean = 0;
        for i = 2:n
            j = B(:,1) == B(i,2);
            if nnz(j) == 1
                % Merge: parent absorbs single child
                B(i,2) = B(j,2);                % Take child's distal node
                B(i,3) = B(j,3);                % Ali Namvar, 2026-02-17 - Use child radius (surviving airway)
                % Recalculate TRUE Euclidean length
                B(i,7) = sqrt(sum((N(B(i,1),2:4) - N(B(i,2),2:4)).^2));
                B(j,:) = [];                    % Remove child branch
                flg_clean = 1;
                break
            end
        end
    end
end

if ~isempty(early_term_log)
    fprintf('Total spurious terminals removed: %d\n', size(early_term_log,1));
else
    fprintf('No spurious terminals found.\n');
end
fprintf('===================================\n\n');

%% ============================================================
%% FINAL CLEANUP AND REINDEXING
%% ============================================================

% Final recalculation of all tree metrics
nB = size(B,1);
S = calc_so(B(:,1:2)); B(:,5) = S(:,end);
H = calc_ho(B(:,1:2)); B(:,6) = H(:,end);
for i = 1:nB
    B(i,4) = gen_no(B(i,1:2), B(:,1:2));
end

% Recalculate all branch lengths (TRUE Euclidean)
for i = 1:nB
    B(i,7) = sqrt(sum((N(B(i,1),2:4) - N(B(i,2),2:4)).^2));
end

% Remove orphaned nodes and reindex
used_nodes = unique([B(:,1); B(:,2)]);

[~, idx] = ismember(used_nodes, N(:,1));
idx = idx(idx > 0);
N = N(idx, :);

% Make node IDs sequential
old_node_ids = N(:,1);
new_node_ids = (1:size(N,1))';

max_old_id = max(old_node_ids);
old_to_new = zeros(max_old_id, 1);
for k = 1:length(old_node_ids)
    old_to_new(old_node_ids(k)) = new_node_ids(k);
end

N(:,1) = new_node_ids;
B(:,1) = old_to_new(B(:,1));
B(:,2) = old_to_new(B(:,2));

%% Final visual check
hf = figure('Name','Final Visual'); ha = axes(hf);
plot_tree(ha,B,N,B(:,4),L_surfs,B_label{4});
for i = 1:nl
    trisurf(L_surfs{i},'FaceColor','b','FaceAlpha',0.05,'EdgeColor','none');
end
saveas(hf,fullfile(outD,[caseID,'.final_visual.fig']));
close(hf);

%% Save data
lobe_ids = [lobe.val];
save(fullfile(outD,[caseID,'.tree_data.mat']),'B','N','lobe_ids','L_surfs','voxsz');

% Save QC reports
if ~isempty(qc_log)
    qc_table = array2table(qc_log, 'VariableNames', ...
        {'Generation', 'Radius_mm', 'Length_mm', 'Type'});
    writetable(qc_table, fullfile(outD, [caseID, '_standard_QC_report.csv']));
end

if ~isempty(early_term_log)
    early_table = array2table(early_term_log, 'VariableNames', ...
        {'Generation', 'Length_mm', 'Strahler_Diff', 'Sibling_Strahler'});
    writetable(early_table, fullfile(outD, [caseID, '_spurious_terminal_QC_report.csv']));
end

fprintf('QC Summary: Standard=%d, Spurious_terminals=%d\n', ...
    size(qc_log,1), size(early_term_log,1));
fprintf('Final tree: %d branches\n', size(B,1));

warning('on','MATLAB:triangulation:PtsNotInTriWarnId')

%% ============================================================
%% LOCAL HELPER FUNCTIONS
%% ============================================================

function str = bool2str(val)
    if val
        str = 'PASS';
    else
        str = 'FAIL';
    end

function y = gen_no(b, B)
    if b(1,1) == B(1,1)
        y = 1;
    else
        b = B(B(:,2) == b(1,1), :);
        % Ali Namvar, 2026-02-17 - Guard against orphan branch (no parent found)
        if isempty(b)
            warning('gen_no: orphan branch detected (no parent), returning NaN');
            y = NaN;
            return;
        end
        y = 1 + gen_no(b, B);
    end
