% CMI script
function prmVoxelAnalyze(cmiObj)
%   Input: cmiObj = CMIclass object containing current settings
%
% Description:
% Note: working, but assumes you reselect PRM model each time to reset
%  thresh and cutoff image values to the correct numbers. Note: also seems
%  to be leaving values of image indices too high???? FIX
%  
% Revision History:
% 17 January 2014 jlb created

% Check for 4 or more images loaded
if length(cmiObj.img.dims)<4
    disp('ERROR: Need at least 4 images to compare PRMs');
    return
end

% Set PRM options
disp('Please set the PRM options before you start!');
cmiObj.setPRMopts();

% Compute prms using sequential pairs
npair = 0;
for pre=1:2:cmiObj.img.dims(4)
    % Compute PRM
    npair = npair + 1;
    cmiObj.img.calcPRM(pre+1);
    prm(npair,:,:,:) = cmiObj.img.prm.mat;
    % Update PRM to choose next pair
    curthresh = cmiObj.img.prm.thresh;
    curthresh(1:2,1:2) = curthresh(1:2,1:2) + 2;
    cmiObj.img.prm.setOpts('thresh',curthresh);
    curcutoff = cmiObj.img.prm.cutoff;
	curcutoff(1,1:2) = curcutoff(1,1:2) + 2; 
    cmiObj.img.prm.setOpts('cutoff',curcutoff);
end


% Find frequency map - assume only one pair for NOW!!!
prm1 = prm(1,:,:,:);
prm2 = prm(2,:,:,:);
prmF = zeros(cmiObj.img.prm.nprm,cmiObj.img.prm.nprm);
for pp=1:cmiObj.img.prm.nprm
    for qq=1:cmiObj.img.prm.nprm % cycle thru columns of prm
        prmF(pp,qq) = nnz(prm1==qq & prm2==pp);
    end
end

% compute cohen kappa stats
kappa(prmF);

% write out array of label counts for t1 and t2 to .csv file
% writeCSV(fullfile(fdir,sprintf('prm%dKappa.csv',cmiObj.img.prm.nprm)),...
%     prm4Kappa,{'PRMLabelT1','PRMLabelT2','VoxelCount'});
% fprintf('Done writing out prmAllLung4Kappa.csv.\n');