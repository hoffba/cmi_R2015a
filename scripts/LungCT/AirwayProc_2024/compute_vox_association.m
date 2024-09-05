function compute_vox_association( ins, seg, A, treedat )

%%% Script for voxel assocaition 
% edited verison of Yixuan's code by AB for single case testing

%% we take the following inputs from tree_prep / virtual output
% TREEd='R:\CGalban_Lab\LabMembers\AlexanderBell\AirwayProc_2024\Airway Testing\out\1018_20240405\1018_20240405.AirwayProc';
% load([TREEd '\1018_20240405.tree_data.mat']);
% load([TREEd '\Branches.mat']);
% load([TREEd '\BranchNodes.mat']);
% load([TREEd '\Strahler.mat']);
% load([TREEd '\Lobe.mat']);

%% data reformatting for voxel association (shifting/renaming to fit an older function)
segBW = logical(seg & ins>-1000 & ins<-500);                          % Throw away values outside of threshold
[xx,yy,zz] = ind2sub(size(seg),logical(seg));  % Voxel index [ i, j, k ]
insHU = ins(segBW);                                                   % Ins HU value
seg = seg(segBW);                                                     % Lobe segmentation ID
A = A(segBW);                                                         % Input values to associate

BRC = [treedat.B()]




% M = get_M_new(LAB,INS,INS); % [x,y,z,eHU,iHU,lobeID] % faster than get_M
% M = trim_M(M,[-1000 -500]); % remove voxels outside given range for analysis
% f_dHU = M(:,4); % functional values for evaluation, e.g. delHU, PRM

BRC = [Branches(:,2:4) Lobe];
BRO = Nodes(:,2:4);

%% voxel association

% BRC ~ N x 4 matrix of branches: [proxID, distID, radius, lobeID]
% BRO ~ N+1 x 3 matrix of nodes [x,y,z] (mm)
% IMB ~ M x 4 matrix of points from CT [x,y,z,f] (mm), f ~ value to average
% Strahler ~ matrix of Strahler orders for BRC (value 1 == terminal branch)
% r ~ radius for obtaining points (set above)

r=4; % arbitrary, just for testing

IMB = [M(:,1:3) f_dHU f_dHU]; % edited to pass dummy variable, will need fixing for registered data
% Yixuan wrote calcVoxel below w assuption of passing exp and ins data, to
% simueltaneously associate eHU, iHU, here I'm rigging it to do
% only one (iHU) essentially; needs generalizing

disp(['Computing voxel assoc...']);
tic
ASSOC = calcVoxel_AB(BRC,BRO,IMB,Strahler,r,pixDim,B,N);
toc

% Yixuan produces huge number of columns, I believe the assigned values
% from voxels are in columns 7 and 8 here; other columns are based on PRM
% and other metrics using both exp and ins HU, and measuring at different
% depths in the tree I believe

B2=[ASSOC(:,[1 2 7])]; % taking just node IDs and assigned values
% N.B. the 3rd column here is the essenital output, the associated voxel
% values (via mean of nearby voxels to terminal airways, and using mean to
% average up to other airways)

% some further processing to get branches in correct order for use with
% other variables such as Strahler

B3=[zeros(size(B2,1),1) B2]; % for branch ids
for i=1:size(B2,1)
    nods=B2(i,1:2);
    [a,b]=ismember(nods,Branches(:,2:3),"rows");
    B3(i,1)=Branches(b,1); % found branch ID
end
B3=sortrows(B3,1); % sort in order of Branch ID