% airway_processing.m
function airway_processing(fn_seg,fn_airways,procdir)

if nargin<2
    procdir = 'R:\CGalban_Lab\LabMembers\BenHoff\tempDATA\AirwayProc_2024';
    fn_seg = fullfile(procdir,'11001N_INSP_STD_BAY_COPD.njh.lobe_segmentation.nii.gz');
    fn_airways = fullfile(procdir,'11001N_INSP_STD_BAY_COPD_Airways.nii.gz');
end
if (nargin < 3) && ~exist('procdir','var')
    procdir = fileparts(fn_airways);
end

% Case ID
[~,pID] = fileparts(fn_seg);
pID = extractBefore(pID,'.');
disp(['Processing subject: ' pID])

% Create folder for outputs
outD = fullfile(procdir,[pID,'.AirwayProc']);
if ~isfolder(outD)
    mkdir(outD);
end

% Prep tree
p = tree_prep_pipe(pID,fn_seg,fn_airways,outD);

CreateConductingZone(p.B(:,1:3), p.N, p.lobe_surfs,[1 1 1],outD);

