function res = Main_Pipeline_sub(expfname,insfname,procdir)

% Input: catalog = cell array of dicom information in a directory generated
%                   by the "dcmCatalog()" script

% Pipeline consists of:
% 1. 

tt = tic;

%% Check that data exists
if ~exist(expfname,'file')
    warning('File does not exist: %s',expfname);
    return;
end
if ~exist(insfname,'file')
    warning('File does not exist: %s',insfname);
    return;
end
if ~isfolder(procdir)
    mkdir(procdir);
end
fn_log = fullfile(procdir,'pipeline_log.txt');

%% Initialize parameters and results struct
fn_ext = '.nii.gz';
check_EI = false;
img = struct('mat',{[],[]},'info',{[],[]},'label',{[],[]});
% tok = regexp(expfname,'\\([^\\\.]+)\.','tokens');
% res.ID = extractBefore(tok{1}{1},'_');
res.ID = '';
res.Exp_DICOM = '';
res.Ins_DICOM = '';
res.ProcDir = procdir;

%% Load CT images
% Exp
writeLog(fn_log,'Loading EXP image...\n');
[img(1).mat,label,fov,orient,info] = cmi_load(1,[],expfname);
d = size(img(1).mat);
img(1).info = struct('label',label,'fov',fov,'orient',orient,'d',d,'voxsz',fov./d,'natinfo',info);
img(1).info.voxvol = prod(img(1).info.voxsz);
img(1).info.name = img(1).info.label;

% Ins
writeLog(fn_log,'Loading INSP image...\n');
[img(2).mat,label,fov,orient,info] = cmi_load(1,[],insfname);
d = size(img(2).mat);
img(2).info = struct('label',label,'fov',fov,'orient',orient,'d',d,'voxsz',fov./d,'natinfo',info);
img(2).info.voxvol = prod(img(2).info.voxsz);
img(2).info.name = img(2).info.label;

if isfolder(expfname)
    res.Exp_DICOM = expfname;
    res.ID = sprintf('%s_%s',info.meta.PatientID,info.meta.StudyDate);
    check_EI = true;
else
    tok = regexp(expfname,'\\([^\\\.]+)\.','tokens');
    res.ID = tok{1}{1};
    if ismember('_',res.ID)
        res.ID = extractBefore(tok{1}{1},'_');
    end
end
res.ElxDir = fullfile(procdir,sprintf('elastix_%s',res.ID));
if isfolder(insfname)
    res.Ins_DICOM = insfname;
    check_EI = true;
end

%% Fix orientation of image:
svchk = false;
tag = {'Exp','Ins'};
for ii = 1:2
    % Find orientation of shoulder bones to see if permute is needed
    BW = img(ii).mat(:,:,end) > -150;
    %BW = max(img(ii).mat(:,:,round(img(ii).info.d(3)/2):end)>800,[],3);
    prop = regionprops(BW,'Orientation','Area');
    if mod(round(prop([prop.Area]==max([prop.Area])).Orientation/90),2)
        writeLog(fn_log,'Permuting %s\n',tag{ii});
        img(ii).mat = permute(img(ii).mat,[2,1,3]);
        img(ii).info.voxsz = img(ii).info.voxsz([2,1,3]);
        img(ii).info.fov = img(ii).info.fov([2,1,3]);
        img(ii).info.d = img(ii).info.d([2,1,3]);
        img(ii).info.orient = img(ii).info.orient([2,1,3,4],[2,1,3,4]);
        svchk = true;
    end
end
clear BW

%% Identify Exp and Ins using lung volume; used for determining file name
if check_EI
    %   ** Need to double-check in case of mislabel
    
    %% Quick segmentation for Exp/Ins Identification:
    img(1).label = getRespiratoryOrgans(img(1).mat);
    img(2).label = getRespiratoryOrgans(img(2).mat);
    
    if (nnz(img(1).label)*img(1).info.voxvol) > (nnz(img(2).label)*img(2).info.voxvol)
        writeLog(fn_log,'Swapping INS/EXP due to lung volume\n');
        img = img([2,1]);
        tstr = res.Exp_DICOM;
        res.Exp_DICOM = res.Ins_DICOM;
        res.Ins_DICOM = tstr;
    end
    svchk = true;
end
img(1).info.label = [res.ID,'_Exp'];
img(2).info.label = [res.ID,'_Ins'];
img(1).info.name = [res.ID,'_Exp'];
img(2).info.name = [res.ID,'_Ins'];

%% Save images:
if svchk
    %% Save nii.gz files using ID and Tag
    fn_exp = fullfile(procdir,sprintf('%s',res.ID,'.exp',fn_ext));
    saveNIFTI(fn_exp,img(1).mat,img(1).info.label,img(1).info.fov,img(1).info.orient);
    fn_ins = fullfile(procdir,sprintf('%s',res.ID,'.ins',fn_ext));
    saveNIFTI(fn_ins,img(2).mat,img(2).info.label,img(2).info.fov,img(2).info.orient);
end

%% Generate Lung Segmentation [This is when VOI don't exist]
fn_label = cell(2,1);
for itag = 1:2
    fn_label{itag} = fullfile(procdir,sprintf('%s.%s.label%s',res.ID,lower(tag{itag}),fn_ext));
    writeLog(fn_log,'%s segmentation ... ',tag{itag});
    if exist(fn_label{itag},'file')
        writeLog(fn_log,'from file\n');
        img(itag).label = readNIFTI(fn_label{itag});
    else
        writeLog(fn_log,'generating using YACTA\n');
        ydir = fullfile(procdir,['yacta_',img(itag).info.name]);
        if ~isfolder(ydir)
            mkdir(ydir);
        end
        tname = fullfile(ydir,sprintf('%s.mhd',img(itag).info.label));
        cmi_save(0,img(itag).mat,img(itag).info.label,img(itag).info.fov,img(itag).info.orient,tname);
        yacta(tname,'wait');
        tname = dir(sprintf('%s*explabels.mhd',tname));
        img(itag).label = cmi_load(1,[],fullfile(ydir,tname.name));
        % Clean-up
        cln_fnames = [dir(fullfile(ydir,'*.mhd'));dir(fullfile(ydir,'*.raw'));dir(fullfile(ydir,'*.zraw'))];
        for ifn = 1:numel(cln_fnames)
            delete(fullfile(ydir,cln_fnames(ifn).name));
        end
        saveNIFTI(fn_label{itag},img(itag).label,img(itag).info.label,img(itag).info.fov,img(itag).info.orient);
    end
end

%% QC segmentation
writeLog(fn_log,'Generating EXP montage...\n');
ind = 10:10:img(1).info.d(3);
QCmontage('seg',cat(4,img(1).mat(:,:,ind),img(1).label(:,:,ind)),img(1).info.voxsz,...
    fullfile(procdir,sprintf('%s_Montage',img(1).info.label)));
writeLog(fn_log,'Generating INSP montage...\n');
ind = 10:10:img(2).info.d(3);
QCmontage('seg',cat(4,img(2).mat(:,:,ind),img(2).label(:,:,ind)),img(2).info.voxsz,...
    fullfile(procdir,sprintf('%s_Montage',img(2).info.label)));

%% ScatterNet for AT on Exp CT scan
fn_scatnet = fullfile(procdir,sprintf('%s.%s%s',res.ID,'scatnet',fn_ext));
writeLog(fn_log,'Air trapping map ... ');
if exist(fn_scatnet,'file')
    writeLog(fn_log,'from file\n');
    atMap = cmi_load(1,img(1).info.d(1:3),fn_scatnet);
else
    writeLog(fn_log,'generating with ScatNet\n');
    atMap = ScatNet(img(1).mat,logical(img(1).label),0);
    cmi_save(0,atMap,'ScatNet',img(1).info.fov,img(1).info.orient,fn_scatnet);
end

%% Quantify unregistered CT scans
writeLog(fn_log,'Quantifying unregistered statistics\n');
S = CTlung_Unreg('exp',img(1).mat,img(1).info.voxvol,img(1).label,atMap);
clear atMap;
for ilab = 1:numel(S)
    if isempty(S(ilab).tag)
        tstr = '';
    elseif ischar(S(ilab).tag)
        tstr = ['_',S(ilab).tag];
    else
        tstr = sprintf('_%u',S(ilab).tag);
    end
    res.(['Exp_Vol',tstr]) = S(ilab).vol;
    res.(['Exp_HU',tstr]) = S(ilab).mean;
    res.(['Exp_856',tstr]) = S(ilab).exp856;
    res.(['Exp_SNpct',tstr]) = S(ilab).SNpct;
    res.(['Exp_SNmean',tstr]) = S(ilab).SNmean;
end
S = CTlung_Unreg('ins',img(2).mat,img(2).info.voxvol,img(2).label);
for ilab = 1:numel(S)
    if isempty(S(ilab).tag)
        tstr = '';
    elseif ischar(S(ilab).tag)
        tstr = ['_',S(ilab).tag];
    else
        tstr = sprintf('_%u',S(ilab).tag);
    end
    res.(['Ins_Vol',tstr]) = S(ilab).vol;
    res.(['Ins_HU',tstr]) = S(ilab).mean;
    res.(['Ins_950',tstr]) = S(ilab).ins950;
    res.(['Ins_810',tstr]) = S(ilab).ins810;
end

%% Register I2E
fn_reg = fullfile(procdir,sprintf('%s.%s%s',res.ID,'ins.reg',fn_ext));
if exist(fn_reg,'file')
    writeLog(fn_log,'Loading registered INS from file...\n');
    ins_reg = readNIFTI(fn_reg);
else
    writeLog(fn_log,'Performing registration...\n')
    lungreg_BH( img(1).mat, img(1).info, logical(img(1).label),...
                img(2).mat, img(2).info, logical(img(2).label),...
                res.ElxDir, res.ID, false, true);
    % Move registered file out to procdir
    [~,outfn] = fileparts(img(2).info.name);
    outfn = fullfile(res.ElxDir,sprintf('%s_R.nii',outfn));
    % Elastix won't save .nii.gz, so need to resave
    ins_reg = readNIFTI(outfn);
    cmi_save(0,ins_reg,'Ins_R',img(1).info.fov,img(1).info.orient,fn_reg);
end
img(2) = [];

%% QC registration
writeLog(fn_log,'Saving Registration Montage ...\n');
ind = 10:10:img(1).info.d(3);
QCmontage('reg',cat(4,ins_reg(:,:,ind),img(1).label(:,:,ind)),img(1).info.voxsz,...
    fullfile(procdir,sprintf('%s_Reg_Montage',res.ID)));

%% PRM calculation
fn_PRM = fullfile(procdir,sprintf('%s.%s%s',res.ID,'prm',fn_ext));
if exist(fn_PRM,'file')
    writeLog(fn_log,'Loading PRM from file...\n');
    prm = int8(readNIFTI(fn_PRM));
else
    writeLog(fn_log,'Calculating PRM...\n');
    [prm,~] = pipeline_PRM(img(1).mat,img(1).info,logical(img(1).label),ins_reg,...
        fullfile(procdir,sprintf('%s_PRM_Scatter',res.ID)));

    % Save PRM
    writeLog(fn_log,'Saving PRM as NIFTI ... ');
    stat = cmi_save(0,prm,{'PRM'},img.info.fov,img.info.orient,fn_PRM);
    if stat
        writeLog(fn_log,'  PRM saved\n');
    else
        writeLog(fn_log,'  Could not save PRM to file.\n');
    end
end
clear ins_reg;

%% map full PRM values (1:10) to (norm,fsad,emph,pd,ns)
prmlabel = {'Norm', 'fSAD', 'Emph', 'PD',       'NS';...
            [1,2],  3,      [4,5],  [8,9,10],   6    };

%% QC PRM
writeLog(fn_log,'Generating PRM Montage ...\n');
QCmontage('prm',cat(4,img.mat(:,:,ind),double(prm(:,:,ind))),...
    img.info.voxsz,fullfile(procdir,sprintf('%s_PRM_Montage',res.ID)));
    
%% Tabulate PRM results:
writeLog(fn_log,'Tabulating PRM results...\n');
BW = logical(img.label);
ulab = unique(img.label(BW));
nlab = numel(ulab);
tstr = '';
for ilab = 1:(nlab+(nlab>1)) % Loop over Whole-lung then R/L
    % Grab sub-region
    if ilab>1
        BW = img.label == (ulab(ilab-1));
        tstr = sprintf('_%u',ulab(ilab-1));
    end
    np = nnz(BW); % Normalize by number in mask
    % 10-color PRM:
    for iprm = 1:10
        res.([sprintf('PRM_%u',iprm),tstr]) = nnz(prm(BW)==iprm)/np*100;
    end
    % 4 color PRM:
    for iprm = 1:size(prmlabel,2)
        res.(['PRM_',prmlabel{1,iprm},tstr]) = nnz(ismember(prm(BW),prmlabel{2,iprm}))/np*100;
    end
end
    
%% Calculate tPRM
prmlabel = ["norm","fsad","emph","pd"];
mflabel = ["v","s","b","x"];
fn_tprm = fullfile(procdir,...
    string(res.ID)+".tprm."+prmlabel+"."+mflabel'+string(fn_ext));
if all(cellfun(@(x)exist(x,'file'),fn_tprm))
    writeLog(fn_log,'Loading tPRM from files ...\n');
    clear prm;
    for iprm = 1:numel(prmlabel)
        for imf = 1:numel(mflabel)
            writeLog(fn_log,'   %s - %s\n',prmlabel(iprm),mflabel(imf))
            tprm = readNIFTI(fn_tprm(imf,iprm));
            for ilab = 1:(nlab+(nlab>1))
                if ilab == 1
                    BW = logical(img.label);
                    tstr = '';
                else
                    BW = img.label == ulab(ilab-1);
                    tstr = sprintf('_%u',ulab(ilab-1));
                end
                res.(strcat('tPRM_',prmlabel(iprm),upper(mflabel(imf)),tstr)) = mean(tprm(BW));
            end
        end
    end
else
    t = tic;
    writeLog(fn_log,'Generating tPRM ...\n');

    %% Calculate MF values
    p = minkowskiFun(prm,'thresh',1:4,...
                         'tmode','==',...
                         'n',10*ones(1,3),...
                         'gridsp',5,...
                         'voxsz',img.info.voxsz,...
                         'mask',logical(img.label),...
                         'prog',0);
    clear prm;

    %% Interpolate to maps
    for ithresh = 1:size(p.MF,1)
        for imf = 1:size(p.MF,2)
            writeLog(fn_log,'   %s - %s\n',prmlabel(ithresh),mflabel(imf));
            
            % Interpolate to image space
            tstr = prmlabel(ithresh) + '.' + mflabel(imf);
            writeLog(fn_log,'       Interpolating\n');
            tprm = grid2img(p.MF(ithresh,imf,:),p.ind,p.mask,3,1);

            % Save tPRM image
            writeLog(fn_log,'       Saving NIFTI\n');
            cmi_save(0,tprm,{char(tstr)},img.info.fov,img.info.orient,char(fn_tprm(imf,ithresh)));

            % Tabulate statistics
            writeLog(fn_log,'       Tabulating means\n');
            for ilab = 1:(nlab+(nlab>1))
                if ilab==1
                    BW = logical(img.label);
                    fstr = '';
                else
                    BW = img.label==ulab(ilab-1);
                    fstr = sprintf('_%u',ulab(ilab-1));
                end
                res.('tPRM_'+prmlabel(ithresh)+'_'+upper(mflabel(imf))+fstr) = mean(tprm(BW));
            end
        end
    end
    writeLog(fn_log,'... tPRM complete (%s)\n',datestr(duration(0,0,toc(t)),'HH:MM:SS'));

end
    
%% Save Results Table:
writetable(struct2table(res,'AsArray',true),fullfile(procdir,[res.ID,'_Results.csv']));

writeLog(fn_log,'Pipeline total time = %s\n',datestr(duration(0,0,toc(tt)),'HH:MM:SS'))


function writeLog(fn,str,varargin)
% write to command window
fprintf(str,varargin{:});

% write to log file
fid = fopen(fn,'a');
if fid
    fprintf(fid,str,varargin{:});
    fclose(fid);
else
    fprintf('Could not open log file.\n');
end

