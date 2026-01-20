% airway_processing.m
function [B,N,lobe_surf] = airway_processing(pID,seg,A,voxsz,procdir)
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

if ischar(seg)
    seg = cmi_load(1,[],seg);
end
if ischar(A)
    A = cmi_load(1,[],A);
end

%% Create folder for outputs
outD = fullfile(procdir,[pID,'.AirwayProc']);
if ~isfolder(outD)
    mkdir(outD);
end

%% Prep tree
t = tic;
[B,N,lobe_surf] = tree_prep_pipe(pID,seg,A,outD,voxsz);
t = round(toc(t));
fprintf('Tree prep = %d:%02d\n',floor(t/60),rem(t,60));
% B ~ [Node1, Node2, Radius, Generation, Strahler, Horsfield, Length]

%% Simulation
t = tic;
[Bsim,Nsim,~,~,num_real_branches] = CreateConductingZone( B(:,1:3), N, lobe_surf , outD );  % add by Ali namvar - get num_real_branches
t = round(toc(t));
fprintf('Simulation = %d:%02d\n',floor(t/60),rem(t,60));
% Bsim ~ [ID, Node1, Node2, Radius, Gen, Strahler, Horsfield, Lobe]
% num_real_branches = number of real (CT-extracted) branches

%% add by Ali namvar - Tag branches based on real/simulated status using index-based approach
t = tic;
nB = size(Bsim, 1);
tag = zeros(nB, 1);

fprintf('\nTagging branches...\n');
fprintf('  Total branches: %d\n', nB);
fprintf('  Real branches (1:%d): %d\n', num_real_branches, num_real_branches);
fprintf('  Simulated branches (%d:%d): %d\n', num_real_branches+1, nB, nB - num_real_branches);

%% add by Ali namvar - Step 1: Identify real terminal branches
% Real terminals are real branches that have at least one SIMULATED daughter
% (meaning they are the roots of simulated subtrees)
real_terminals = [];
real_nonterminals = [];

for i = 1:num_real_branches
    dist_node = Bsim(i, 3);  % distal node of this real branch
    
    % Find all daughters (branches whose proximal node = this distal node)
    daughters = find(Bsim(:, 2) == dist_node);
    
    if isempty(daughters)
        % No daughters at all - this is a terminal that had no simulation
        % (possibly outside lobe boundaries or other reason)
        real_terminals = [real_terminals; i];
    else
        % Check if any daughters are simulated
        sim_daughters = daughters(daughters > num_real_branches);
        real_daughters = daughters(daughters <= num_real_branches);
        
        if ~isempty(sim_daughters)
            % Has simulated daughters - this is a real terminal
            real_terminals = [real_terminals; i];
        else
            % Only has real daughters - this is a real non-terminal
            real_nonterminals = [real_nonterminals; i];
        end
    end
end

fprintf('  Real terminals (roots of simulated trees): %d\n', length(real_terminals));
fprintf('  Real non-terminals: %d\n', length(real_nonterminals));

%% add by Ali namvar - Step 2: Tag real terminals with negative numbers
for i = 1:length(real_terminals)
    tag(real_terminals(i)) = -i;
end

% Real non-terminals keep tag = 0

%% add by Ali namvar - Step 3: Tag simulated branches with positive numbers
% Each simulated branch gets tagged with the index of its real terminal root

% Build a mapping: for each simulated branch, find which real terminal it belongs to
% We do this by tracing back through parents

for i = (num_real_branches + 1):nB
    % This is a simulated branch
    % Trace back through parents until we find a real branch
    current = i;
    max_iterations = nB;  % Safety limit
    iter = 0;
    
    while current > num_real_branches && iter < max_iterations
        prox_node = Bsim(current, 2);  % proximal node of current branch
        
        % Find parent branch (whose distal node = this proximal node)
        parent = find(Bsim(:, 3) == prox_node, 1);
        
        if isempty(parent)
            % No parent found - orphan branch (shouldn't happen)
            warning('Orphan simulated branch found: %d', i);
            break;
        end
        
        current = parent;
        iter = iter + 1;
    end
    
    if current <= num_real_branches
        % Found the real branch that this simulated branch descends from
        % Find which terminal index this corresponds to
        term_idx = find(real_terminals == current);
        if ~isempty(term_idx)
            tag(i) = term_idx;  % positive = which terminal it belongs to
        else
            % The real branch is not in our terminal list - shouldn't happen
            % but handle gracefully
            warning('Simulated branch %d traces to real branch %d which is not a terminal', i, current);
            tag(i) = 999999;  % Flag as problematic
        end
    else
        % Couldn't trace back to a real branch
        warning('Could not trace simulated branch %d to real root', i);
        tag(i) = 999998;  % Flag as problematic
    end
end

%% Validate tagging
n_real_nonterminal = sum(tag == 0);
n_real_terminal = sum(tag < 0);
n_simulated = sum(tag > 0);
n_problematic = sum(tag >= 999998);

fprintf('\nTag validation:\n');
fprintf('  Tag = 0 (real non-terminal): %d\n', n_real_nonterminal);
fprintf('  Tag < 0 (real terminal): %d\n', n_real_terminal);
fprintf('  Tag > 0 (simulated): %d\n', n_simulated);

if n_problematic > 0
    warning('%d branches have problematic tags (could not trace to root)', n_problematic);
end

% Sanity checks
expected_real_nonterminal = length(real_nonterminals);
expected_real_terminal = length(real_terminals);
expected_simulated = nB - num_real_branches;

if n_real_nonterminal ~= expected_real_nonterminal
    warning('Real non-terminal count mismatch: got %d, expected %d', ...
        n_real_nonterminal, expected_real_nonterminal);
end

if n_real_terminal ~= expected_real_terminal
    warning('Real terminal count mismatch: got %d, expected %d', ...
        n_real_terminal, expected_real_terminal);
end

if n_simulated ~= expected_simulated
    warning('Simulated count mismatch: got %d, expected %d', ...
        n_simulated, expected_simulated);
end

t = round(toc(t));
fprintf('Branch tagging = %d:%02d\n',floor(t/60),rem(t,60));

%% Assemble output
B = [Bsim, tag];
N = Nsim;

%% Save tree
B_label = {'ID','N_Prox','N_Dist','Radius','Gen','Strahler','Horsfield','Lobe','Tag'};
dims = size(seg);

%% add by Ali namvar - Save tag statistics
tag_stats.total_branches = nB;
tag_stats.num_real_branches = num_real_branches;
tag_stats.real_nonterminal = n_real_nonterminal;
tag_stats.real_terminal = n_real_terminal;
tag_stats.simulated = n_simulated;
tag_stats.real_terminal_indices = real_terminals;

save(fullfile(outD,[pID,'_AirwayTreeSim.mat']),'dims','voxsz','B','N','B_label','lobe_surf','tag_stats','num_real_branches');

fprintf('\nSaved airway tree to: %s\n', fullfile(outD,[pID,'_AirwayTreeSim.mat']));
fprintf('Final tree: %d branches (%d real, %d simulated)\n', nB, num_real_branches, nB-num_real_branches);
