function Main_Pipeline(catalog)

% Input: catalog = cell array of dicom information in a directory generated
%                   by the "dcmCatalog()" script

% Pipeline consists of:
% 1. Identify image key information and save as .nii
% 2. Perform elastix registration

%% Determine data for processing:
% cases is structure array:
%   cases(i).PatientID
%           .StudyDate
%           .Scans - Structure containing scan info for selected data
if nargin==0
    catalog = [];
end
[cases, home_pwd] = catalog_select_2(catalog);

%% prepare Data Table with Dummy values
dataTable = table;

% Find RegClass object in base workspace
regObj = getObjectsFromBase('RegClass',1);
% If non found, create a temporary one without GUI
if isempty(regObj)
    regObj = RegClass(0);
end

%% Loop over cases
N = length(cases);
fn_ext = '.nii.gz';
h1 = waitbar(0, 'Analyze Exp and Ins data');
for i = 1:N
    ID = cases(i).PatientName;
    waitbar(i/N,h1,[num2str(round(100*(i-1)/N,1)),'% Complete: Load Data and Masks for ',ID])
        
        %% Establish relevant filenmes:
        procdir = fullfile(home_pwd,'ProcessedData',ID);
        fname_Exp = sprintf('%s_Exp',ID);
        fname_Ins = sprintf('%s_Ins',ID);
        fname_Exp_Label = sprintf('%s_Exp_Label',ID);
        fname_Ins_Label = sprintf('%s_Ins_Label',ID);
        fname_ScatNet = sprintf('%s_SNmap',ID);
        
        %% Initialize data structure
%         data = struct('tag',{'Exp','Ins'},'mat',{[],[]},'voi',{[],[]},'info',{[],[]});
        
        %% Read selected DICOM data:
        fprintf('\nReading image data from file ... ID = %s\n',ID)
        if ~(exist(fullfile(procdir,[fname_Exp,fn_ext]),'file') && exist(fullfile(procdir,[fname_Ins,fn_ext]),'file'))
            fprintf('   ... from DICOM\n');
            data = Step01_dcm2nii(procdir,{cases(i).Scans(:).Directory},{cases(i).Scans(:).Tag});
            check_EI = true;
        else
            fprintf('   ... from NiFTi\n');
            
            data(1).img.info = niftiinfo(fullfile(procdir,[fname_Exp,fn_ext]));
            data(1).img.mat = niftiread(data(1).info); 
            data(1).tag = 'Exp';
            
            data(2).img.info = niftiinfo(fullfile(procdir,[fname_Ins,fn_ext])); 
            data(2).img.mat = niftiread(data(2).info); 
            data(2).tag = 'Ins';
            
            check_EI = false;
        end
        
        %% Generate Lung Segmentation [This is when VOI don't exist]
        if ~(exist(fullfile(procdir,fname_Exp_Label),'file') && exist(fullfile(procdir,fname_Ins_Label),'file'))
            fprintf('   Generating VOIs\n')
            [data(1).voi.mat, data(1).voi.info] = Step02_segLungHuman_cjg(1,data(1).img.mat,data(1).img.info, procdir);
            [data(2).voi.mat, data(2).voi.info] = Step02_segLungHuman_cjg(1,data(2).img.mat,data(2).img.info, procdir);
        else
            fprintf('   Reading VOIs from file\n')
            data(1).voi.info = niftiinfo(fullfile(procdir,fname_Exp_Label));
            data(1).voi.mat = niftiread(data(1).voi.info);
            data(2).voi.info = niftiinfo(fullfile(procdir,fname_Ins_Label));
            data(2).voi.mat = niftiread(data(2).voi.info);
        end

        %% Identify Exp and Ins using lung volume; used for determining file name
        if check_EI
            %   ** Need to double-check in case of mislabel
            idx_EIV = Step03_Identify_EI(data);
            data = data(idx_EIV);
        end
        
        %% Save nii.gz files using ID and Tag
        Step04_Rename_EI(procdir,fname_Exp,fname_Exp_Label,fname_Ins,fname_Ins_Label,data);
        
        %% Quantify unregistered CT scans
        fprintf('\n   Quantifying unregistered statistics\n');
        [expData,insData] = Step05_UnRegLungAnalysis(procdir, fname_ScatNet, data);

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
        lungreg_BH( fullfile(procdir,[fname_Exp,fn_ext]), fullfile(procdir,[fname_Ins,fn_ext]),...
            fullfile(procdir,[fname_Exp_Label,fn_ext]), fullfile(procdir,[fname_Ins_Label,fn_ext]), regObj);
end
delete(h1);

%% Assign data to base workspace
assignin('base', 'dataTable', dataTable)