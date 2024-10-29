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

% Create folder for outputs
outD = fullfile(procdir,[pID,'.AirwayProc']);
if ~isfolder(outD)
    mkdir(outD);
end

% Prep tree
[B,N,lobe_surf] = tree_prep_pipe(pID,seg,A,outD,voxsz);

% Simulation
[Bsim,Nsim,~,~] = CreateConductingZone( B(:,1:3), N, lobe_surf , outD );

% Tag real vs simulated branches
Nreal = Nsim(ismember(Nsim(:,2:4),N(:,2:4),'rows'),1); % Find original nodes
B = [ Bsim , sum(ismember(Bsim(:,2:3),Nreal),2)==2 ]; % Find branches between original nodes
N = Nsim;

% Save tree
B_label = {'ID','N_Prox','N_Dist','Radius','Gen','Strahler','Horsfield','Lobe','Real','Terminal'};
dims = size(seg);
save(fullfile(outD,[pID,'_AirwayTreeSim.mat']),'dims','voxsz','B','N','B_label','lobe_surf');












return
%% Compute voxel associations with branches
segBW = seg & ins>-1000 & ins<-500;         % Throw away values outside of threshold
[xi,yi,zi] = ind2sub(dim,find(segBW));      % Voxel index [ i, j, k ]
nv = size(img,4);
V = reshape(img(repmat(segBW,1,1,1,nv)),[],nv);

% r ~ radius for obtaining points (set above)
r=4; % arbitrary, just for testing

%% temp set up for old version of calVoxel
t = tic;
ASSOC = calcVoxel_AB( Bsim(:,[2,3,4,8]) , Nsim(:,2:4) , [xi,yi,zi,V,V] , Bsim(:,5) , r , voxsz , B , N );
toc(t);

% Yixuan produces huge number of columns, I believe the assigned values
% from voxels are in columns 7 and 8 here; other columns are based on PRM
% and other metrics using both exp and ins HU, and measuring at different
% depths in the tree I believe

nB = size(B,1);
B = [zeros(nB,1), ASSOC(:,[1 2 7])]; % taking just node IDs and assigned values
% N.B. the 3rd column here is the essenital output, the associated voxel
% values (via mean of nearby voxels to terminal airways, and using mean to
% average up to other airways)

% some further processing to get branches in correct order for use with
% other variables such as Strahler

% B3 = [zeros(size(B2,1),1) B2]; % for branch ids
for i = 1:nB
    [~,b] = ismember(B(i,2:3),B{:,2:3},"rows");
    B(i,1) = B.ID(b); % found branch ID
end
B = sortrows(B,1); % sort in order of Branch ID

% Save to procdir
save(fullfile(outD,'ASSOC.mat'),'B');
