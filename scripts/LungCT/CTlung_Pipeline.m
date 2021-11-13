function results = CTlung_Pipeline(info,basedir,measures)
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

%% Initialize processing structure:
fn_base = sprintf('%s_%s',info.PatientID,info.StudyDate);
ct.exp = ImageClass;
ct.ins = ImageClass;
ct.procdir = fullfile(basedir,'ProcessedData',fn_base);
ct.elxdir = fullfile(ct.procdir,sprintf('elxreg_%s',fn_base));
ct.ext = '.nii.gz';

%% Check dependencies and previously saved data
% opts = CTlung_CheckDependencies(ct,measures);

%% Read selected DICOM data:
check_EI = CTlung_ReadImages(ct);

%% Generate Lung Segmentation [This is when VOI don't exist]
CTlung_Segmentation();
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
