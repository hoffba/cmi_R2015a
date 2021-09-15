function Main_Pipeline_Batch(catalog)

if nargin==0
    catalog = [];
elseif ~(ischar(catalog) && exist(catalog,'file'))
    error('Invalid input data catalog');
end

[cases, home_pwd] = catalog_select_2(catalog);

%% Start batch process for each case:
for i = 1:length(cases)
    
end

%% Wait for completion, then compile results for table
for i = 1:length(cases)
end

end


function results = Batch_Pipeline(ID,procdir,dcmnames)

%% Initialize processing structure:
data = struct(...
    'id',{ID},...
    'img',{[ImageClass;ImageClass]},...
    'dcmnames',{dcmnames},...
    'procdir',fullfile(procdir,ID),...
    'elxdir',fullfile(procdir,ID,sprintf('elxreg_%s',fname_Ins)),...
    'fname',struct('exp',sprintf('%s_Exp',ID),...
                   'ins',sprintf('%s_Ins',ID),...
                   'exp_label',sprintf('%s_Exp_Label',ID),...
                   'ins_label',sprintf('%s_Ins_Label',ID),...
                   'scatnet',sprintf('%s_SNmap',ID)),...
    'ext','.nii.gz',...
    'results',{});

%% Read images
CTlung_ReadImages(data,dcmnames);
        
%% Determine lung segmenation
CTlung_Segmentation(data,3)

%% Identify Exp and Ins using lung volume; used for determining file name
if check_EI
    CTlung_checkEI(regObj,fnames,procdir)
end

%% Quantify unregistered CT scans


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

end
