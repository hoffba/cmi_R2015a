function res = Main_Pipeline_sub(expfname,insfname,procdir)

% Input: catalog = cell array of dicom information in a directory generated
%                   by the "dcmCatalog()" script

% Pipeline consists of:
% 1. 

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

%% Initialize parameters and results struct
fn_ext = '.nii.gz';
check_EI = false;
img = struct('mat',{[],[]},'info',{[],[]},'label',{[],[]});
tok = regexp(expfname,'\\([^\\\.]+)\.','tokens');
res.ID = extractBefore(tok{1}{1},'_');
res.ID = '';
res.Exp_DICOM = '';
res.Ins_DICOM = '';
res.ProcDir = procdir;

%% Load CT images
% Exp
fprintf('Loading EXP image...\n');
[img(1).mat,label,fov,orient,info] = cmi_load(1,[],expfname);
d = size(img(1).mat);
img(1).info = struct('label',label,'fov',fov,'orient',orient,'d',d,'voxsz',fov./d,'natinfo',info);
img(1).info.voxvol = prod(img(1).info.voxsz);
img(1).info.name = img(1).info.label;

% Ins
fprintf('Loading INSP image...\n');
[img(2).mat,label,fov,orient,info] = cmi_load(1,[],insfname);
d = size(img(2).mat);
img(2).info = struct('label',label,'fov',fov,'orient',orient,'d',d,'voxsz',fov./d,'natinfo',info);
img(2).info.voxvol = prod(img(2).info.voxsz);
img(2).info.name = img(2).info.label;

if isfolder(expfname)
    res.Exp_DICOM = expfname;
    res.ID = info.PatientName;
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
    res.Exp_DICOM = insfname;
    check_EI = true;
end

%% Fix orientation of image:
svchk = false;
tag = {'Exp','Ins'};
for ii = 1:2
    % Find orientation of shoulder bones to see if permute is needed
    BW = max(img(ii).mat(:,:,round(img(ii).info.d(3)/2):end)>800,[],3);
    prop = regionprops(BW,'Orientation','Area');
    if mod(round(prop([prop.Area]==max([prop.Area])).Orientation/90),2)
        fprintf('Permuting %s\n',tag{ii});
        img(ii).mat = permute(img(ii).mat,[2,1,3]);
        img(ii).info.voxsz = img(ii).info.voxsz([2,1,3]);
        img(ii).info.fov = img(ii).info.fov([2,1,3]);
        img(ii).info.d = img(ii).info.d([2,1,3]);
        img(ii).info.orient = img(ii).info.orient([2,1,3,4],[2,1,3,4]);
        svchk = true;
    end
end

%% Quick segmentation for Exp/Ins Identification:
img(1).label = getRespiratoryOrgans(img(1).mat);
img(2).label = getRespiratoryOrgans(img(2).mat);

%% Identify Exp and Ins using lung volume; used for determining file name
if check_EI
    %   ** Need to double-check in case of mislabel
    if (nnz(img(1).label)*img(1).info.voxvol) > (nnz(img(2).label)*img(2).info.voxvol)
        fprintf('Swapping INS/EXP due to lung volume\n');
        img = img([2,1]);
        tstr = res.Exp_DICOM;
        res.Exp_DICOM = res.Ins_DICOM;
        res.Ins_DICOM = tstr;
    end
    svchk = true;
end

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
    fprintf('%s segmentation ... ',tag{itag});
    if exist(fn_label{itag},'file')
        fprintf('from file\n');
        img(itag).label = readNIFTI(fn_label{itag});
    else
        fprintf('generating using YACTA\n');
        ydir = fullfile(procdir,['yacta_',img(itag).info.label]);
        if ~isfolder(ydir)
            mkdir(ydir);
        end
        tname = fullfile(ydir,sprintf('%s.mhd',img(itag).info.label));
        cmi_save(0,img(itag).mat,img(itag).info.label,img(itag).info.fov,img(itag).info.orient,tname);
        yacta(tname,'wait');
        tname = dir(sprintf('%s*.mhd',tname));
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
% fprintf('Generating EXP montage...\n');
% ind = 10:10:img(1).info.d(3);
% QCmontage('seg',cat(4,img(1).mat(:,:,ind),img(1).label(:,:,ind)),img(1).info.voxsz,...
%     fullfile(procdir,sprintf('%s.Montage.tiff',img(1).info.label)));
% fprintf('Generating INSP montage...\n');
% ind = 10:10:img(2).info.d(3);
% QCmontage('seg',cat(4,img(2).mat(:,:,ind),img(2).label(:,:,ind)),img(2).info.voxsz,...
%     fullfile(procdir,sprintf('%s.Montage.tiff',img(2).info.label)));

%% ScatterNet for AT on Exp CT scan
fn_scatnet = fullfile(procdir,sprintf('%s.%s%s',res.ID,'scatnet',fn_ext));
fprintf('Air trapping map ... ');
if exist(fn_scatnet,'file')
    fprintf('from file\n');
    atMap = cmi_load(1,img(1).info.d(1:3),fn_scatnet);
else
    fprintf('generating with ScatNet\n');
    atMap = ScatNet(img(1).mat,logical(img(1).label),0);
    cmi_save(0,atMap,'ScatNet',img(1).info.fov,img(1).info.orient,fn_scatnet);
end

%% Quantify unregistered CT scans
fprintf('Quantifying unregistered statistics\n');
S = CTlung_Unreg('exp',img(1).mat,img(1).info.voxvol,img(1).label,atMap);
clear atMap;
for ilab = 1:numel(S)
    if isempty(S(ilab).tag)
        tstr = '';
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
    fprintf('Loading registered INS from file...\n');
    ins_reg = readNIFTI(fn_reg);
else
    fprintf('Performing registration...\n')
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
fprintf('Saving Registration Montage ...\n');
ind = 10:10:img(1).info.d(3);
QCmontage('reg',cat(4,ins_reg(:,:,ind),img(1).label(:,:,ind)),img(1).info.voxsz,...
    fullfile(procdir,sprintf('%s_Reg_Montage',res.ID)));

%% PRM calculation
fn_PRM = fullfile(procdir,sprintf('%s.%s%s',res.ID,'prm',fn_ext));
if exist(fn_PRM,'file')
    fprintf('Loading PRM from file...\n');
    prm = int8(readNIFTI(fn_PRM));
else
    fprintf('Calculating PRM...\n');
    [prm,~] = pipeline_PRM(img(1).mat,img(1).info,logical(img(1).label),ins_reg,...
        fullfile(procdir,sprintf('%s_PRM_Scatter',res.ID)));

    % Save PRM
    fprintf('Saving PRM as NIFTI ... ');
    stat = cmi_save(0,prm,{'PRM'},img.info.fov,img.info.orient,fn_PRM);
    if stat
        fprintf('  PRM saved\n');
    else
        fprintf('  Could not save PRM to file.\n');
    end
end
clear ins_reg;
% map full PRM values (1:10) to (norm,fsad,emph,pd,ns)
prmlabel = {'Norm','fSAD','Emph','PD','NS'};
prm(ismember(prm,1:2)) = 1; % Norm
prm(prm==3)            = 2; % fSAD
prm(ismember(prm,4:5)) = 3; % Emph
prm(prm==6)            = 5; % NS
prm(ismember(prm,8:10))= 4; % PD

% QC PRM
fprintf('Generating PRM Montage ...\n');
QCmontage('prm',cat(4,img.mat(:,:,ind),double(prm(:,:,ind))),...
    img.info.voxsz,fullfile(procdir,sprintf('%s_PRM_Montage',res.ID)));
    
%% Tabulate PRM results:
fprintf('Tabulating PRM results...\n');
BW = logical(img.label);
ulab = unique(img.label(BW));
nlab = numel(ulab);
tstr = '';
for ilab = 1:(nlab+(nlab>1))
    if ilab>1
        BW = img.label == (ulab(ilab-1));
        tstr = sprintf('_%u',ulab(ilab-1));
    end
    np = nnz(BW); % Normalize by number in mask or nnz in PRM?????????
    for iprm = 1:numel(prmlabel)
        res.(['PRM_',prmlabel{iprm},tstr]) = nnz(prm(BW)==iprm)/np*100;
    end
end
    
%% Calculate tPRM
prmlabel = ["norm","fsad","emph","pd"];
mflabel = ["v","s","b","x"];
fn_tprm = fullfile(procdir,...
    string(res.ID)+".tprm."+prmlabel+"."+mflabel'+string(fn_ext));
if all(cellfun(@(x)exist(x,'file'),fn_tprm))
    fprintf('Loading tPRM from files ...\n');
    clear prm;
    for iprm = 1:numel(prmlabel)
        for imf = 1:numel(mflabel)
            fprintf('   %s - %s\n',prmlabel(iprm),mflabel(imf))
            tprm = readNIFTI(fn_tprm(imf,iprm));
            for ilab = 1:(nlab+(nlab>1))
                if ilab == 1
                    BW = logical(img.label);
                    tstr = '';
                else
                    BW = img.label == ulab(ilab-1);
                    tstr = sprintf('_%u',ulab(ilab-1));
                end
                res.(['tPRM_',prmlabel(iprm),upper(mflabel(imf)),tstr]) = mean(tprm(BW));
            end
        end
    end
else
    t = tic;
    fprintf('Generating tPRM ...\n');

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
            fprintf('   %s - %s\n',prmlabel(ithresh),mflabel(imf));
            
            % Interpolate to image space
            tstr = prmlabel(ithresh) + '.' + mflabel(imf);
            fprintf('       Interpolating\n');
            tprm = grid2img(p.MF(ithresh,imf,:),p.ind,p.mask,3,1);

            % Save tPRM image
            fprintf('       Saving NIFTI\n');
            cmi_save(0,tprm,{char(tstr)},img.info.fov,img.info.orient,char(fn_tprm(imf,ithresh)));

            % Tabulate statistics
            fprintf('       Tabulating means\n');
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
    fprintf('... tPRM complete (%s)\n',datestr(duration(0,0,toc(t)),'HH:MM:SS'));

end
    
%% Save Results Table:
writetable(struct2table(res,'AsArray',true),fullfile(procdir,[res.ID,'_Results.csv']));

