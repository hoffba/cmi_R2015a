% airway_processing.m
function [B,N,lobe_surf] = airway_processing(pID,seg,A,voxsz,procdir)
%% function airway_processing(pID,ins,seg,A,voxsz,img,fcn,procdir)
%   Inputs: pID =       string ID of subject/case (base filename)
%           ins =       inspiratory image
%           seg =       inspiratory lung segmentation
%           A =         airways segmentation
%           voxsz =     voxel dimensions
%           img =       4D image matrix for association with airway tree
%           label =     cell array of labels for each input img
%           fcn =       cell array with function handles to apply to association window
%                           ** functions must have single output and vector input
%           procdir =   directory to save results to
%   Outputs: B = [ID, N_prox, N_dist, Radius, Gen, Strahl, Hors, Lobe, Real, Terminal]
%            N = [ID, x, y, z]

%% Create folder for outputs
outD = fullfile(procdir,[pID,'.AirwayProc']);
if ~isfolder(outD)
    mkdir(outD);
end

%% Prep tree
t = tic;
[B,N,lobe_surf] = tree_prep_pipe(pID,seg,A,outD,voxsz);
t = round(toc(t));
fprintf('Tree prep = %d:%d\n',floor(t/60),rem(t,60));
% B ~ [Node1, Node2, Radius, Horsfield, Strahler, Generation]

%% Simulation
t = tic;
[Bsim,Nsim,~,~] = CreateConductingZone( B(:,1:3), N, lobe_surf , outD );
t = round(toc(t));
fprintf('Simulation = %d:%d\n',floor(t/60),rem(t,60));
% Bsim ~ [ID, Node1, Node2, Radius, ...]

%% Tag simulated branches
%  0    = Real branch
% -1:-N = Real terminal
%  1:N  = Simulated (# associated with -# for real terminal)
t = tic;
tag = zeros(size(Bsim,1),1);

% Find real terminal branches
Nr_term = N( ismember(N(:,1),B(B(:,5)==1,2)) , 2:4 );   % xyz
Ns_term = Nsim(ismember(Nsim(:,2:4),Nr_term,'rows'),1); % ID
Br_term = find(ismember( Bsim(:,3) , Ns_term ));
nterm = numel(Br_term);
tag(Br_term) = -(1:nterm);

% Label each simulated branch
for i = 1:nterm
    Bnext = ismember (Bsim(:,2) , Bsim(Br_term(i),3) );
    while nnz(Bnext)
        tag(Bnext) = i;

        % Find daughter branches
        Bnext = ismember( Bsim(:,2) , Bsim(Bnext,3) );
    end
end
B = [ Bsim , tag ];
N = Nsim;
t = round(toc(t));
fprintf('Branch tags = %d:%d\n',floor(t/60),rem(t,60));

%% Save tree
B_label = {'ID','N_Prox','N_Dist','Radius','Gen','Strahler','Horsfield','Lobe','Tag'};
dims = size(seg);
save(fullfile(outD,[pID,'_AirwayTreeSim.mat']),'dims','voxsz','B','N','B_label','lobe_surf');



