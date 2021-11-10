function output = CTlung_Pipeline(ct,measures)
% Inputs: info = [struct] containing case info
%             .ProcDir = [char] base path for data and saving
%             .PatientID = [char]
%             .StudyDate = [char]
%             .Scans = [struct] containing scan info for raw INS/EXP dicoms
%                   ** should have two elements
%         measures = [string array] containing desired outputs
%                   ** defaults to all if empty
%             .
% Output: output = [struct] containing analysis results

% Pipeline consists of:
% 1. 


%% prepare case info
output = [];
procdir = fullfile(ct.BasePath,'ProcessedData',ct.PatientID);
fn_base = sprintf('%s_%s',ct.PatientID,ct.StudyDate);

%% Initialize processing structure:
ct.exp = ImageClass;
ct.ins = ImageClass;
ct.elxdir = fullfile(ct.ProcDir,sprintf('elxreg_%s',fn_base));
% info.fn = struct('exp',sprintf('%s_Exp',fn_base),...
%                  'ins',sprintf('%s_Ins',fn_base),...
%                  'exp_label',sprintf('%s_Exp_Label',fn_base),...
%                  'ins_label',sprintf('%s_Ins_Label',fn_base),...
%                  'scatnet',sprintf('%s_SNmap',fn_base),...
%                  'csa',sprintf('%s_CSA_skel.mat',fn_base),...
%                  'prm',sprintf('%s_PRM',fn_base),...
%                  'tprm',sprintf('%s_',fn_base));
ct.ext = '.nii.gz';
ct.Results = [];

%% Check dependencies and previously saved data
opts = CTlung_CheckDependencies(ct,measures);

%% Read selected DICOM data:
fprintf('\nReading image data from file ... ID = %s\n',ID)
stat = CTlung_ReadImages(ct,opts);

if ~(exist(fullfile(procdir,[fname_Exp,fn_ext]),'file') && exist(fullfile(procdir,[fname_Ins,fn_ext]),'file'))
    fprintf('   ... from DICOM\n');
    regObj.cmiObj(1).loadImg(0,cases(i).Scans(strcmp({cases(i).Scans(:).Tag},'Exp')).Directory,procdir,fname_Exp);
    regObj.cmiObj(2).loadImg(0,cases(i).Scans(strcmp({cases(i).Scans(:).Tag},'Ins')).Directory,procdir,fname_Ins);
    check_EI = true;
else
    fprintf('   ... from NiFTi\n');
    fprintf('   ... Reading Exp\n');
    regObj.cmiObj(1).loadImg(0,fullfile(procdir,[fname_Exp,fn_ext]),procdir,fname_Exp);
    fprintf('   ... Reading Ins\n');
    regObj.cmiObj(2).loadImg(0,fullfile(procdir,[fname_Ins,fn_ext]),procdir,fname_Ins);
    check_EI = false;
end

%% Generate Lung Segmentation [This is when VOI don't exist]
if ~(exist(fullfile(procdir,[fname_Exp_Label,fn_ext]),'file') && exist(fullfile(procdir,[fname_Ins_Label,fn_ext]),'file'))
    fprintf('   Generating VOIs\n')
    tmask = Step02_segLungHuman_cjg(1,regObj.cmiObj(1).img.mat,fname_Exp_Label, procdir);
    regObj.cmiObj(1).img.mask.merge('replace',tmask);
    tmask = Step02_segLungHuman_cjg(1,regObj.cmiObj(2).img.mat,fname_Ins_Label, procdir);
    regObj.cmiObj(2).img.mask.merge('replace',tmask);
else
    fprintf('   Reading VOIs from file\n')
    regObj.cmiObj(1).loadMask(fullfile(procdir,[fname_Exp_Label,fn_ext]));
    regObj.cmiObj(2).loadMask(fullfile(procdir,[fname_Ins_Label,fn_ext]));
end

%% Identify Exp and Ins using lung volume; used for determining file name
if check_EI
    %   ** Need to double-check in case of mislabel
    if nnz(regObj.cmiObj(1).img.mask.mat) > nnz(regObj.cmiObj(2).img.mask.mat)
        regObj.swapCMIdata;
    end
    %% Save nii.gz files using ID and Tag
    regObj.cmiObj(1).img.saveImg(1,fullfile(procdir,[fname_Exp,fn_ext]),1);
    regObj.cmiObj(2).img.saveImg(1,fullfile(procdir,[fname_Ins,fn_ext]),1);
    regObj.cmiObj(1).img.saveMask(fullfile(procdir,[fname_Exp_Label,fn_ext]));
    regObj.cmiObj(2).img.saveMask(fullfile(procdir,[fname_Ins_Label,fn_ext]));
end


%% Quantify unregiste44444red CT scans
fprintf('\n   Quantifying unregistered statistics\n');
[expData,insData] = Step05_UnRegLungAnalysis(procdir, fname_ScatNet, regObj);

%% Save data to Table
dataTable.ID(i,:) = ID;

dataTable.Exp_dFile(i,:) = {fname_Exp};
dataTable.Ins_dFile(i,:) = {fname_Ins};

dataTable.Exp_mFile(i,:) = {fname_Exp_Label};
dataTable.Ins_mFile(i,:) = {fname_Ins_Label};

dataTable.ExpVol(i,1) = expData(1); dataTable.ExpHU(i,1) = expData(2); dataTable.Exp856(i,1) = expData(3);
dataTable.SNperc(i,1) = expData(4); dataTable.SNmean(i,1) = expData(5);

dataTable.InsVol(i,1) = insData(1); dataTable.InsHU(i,1) = insData(2); dataTable.Ins950(i,1) = insData(3);
dataTable.Ins810(i,1) = insData(4);

%% Register I2E
lungreg_BH(ID,elxdir,regObj);
