function Main_Pipeline_GL(varargin)
% Start CT lung pipeline script (Main_Pipeline_sub) with GUI selection for
% cases based on DCMcatalog.csv
%
% Main_Pipeline_select( dcmpath, svpath )
%       * runs cases using full reg and default local cluster settings
% Main_Pipeline_select( dcmpath, svpath, 'quickreg' )
%       * skips last step of registration for quick results
% Main_Pipeline_select( dcmpath, svpath, 'serial' )
%       * runs cases serially in the command line (for debugging)
% Main_Pipeline_select( dcmpath, svpath, 'serialbatch' )
%       * runs cases serially in batch (1-worker)
% Main_Pipeline_select( dcmpath, svpath, 'quickreg', 'serialbatch' )

cluster_profile = {};
flag = false;
opts = struct('dcmpath','',...
              'sv_path','',...
              'quickreg',false);

if nargin
    % Find directory inputs
    pathi = find(cellfun(@(x)ischar(x) && isfolder(x),varargin),2);
    if numel(pathi)
        opts.dcmpath = varargin{pathi(1)};
    end
    if numel(pathi)>1
        opts.sv_path = varargin{pathi(2)};
    end
    if ismember('quickreg',varargin)
        opts.quickreg = true;
    end
end


%% Find files
[cases,opts] = catalog_select_3(opts);

%% Make sure the save directory is on Turbo, for access from GL
[~,TF] = checkTurboPath(opts.savepath);
if TF
    error('Save path must be on Turbo for access from Great Lakes.');
end

%% Determine batch setup
if ismember('serial',varargin)
    flag = true;
elseif ismember('serialbatch',varargin)
    cluster_profile = {'Profile','cmi_1'};
    if ~ismember('cmi_1',parallel.clusterProfiles)
        mls_fname = fullfile(fileparts(which('cmi')),'cmi_1.mlsettings');
        if exist(mls_fname,'file')
            fprintf('Importing cluster profile: %s\n',mls_fname);
            parallel.importProfile(mls_fname);
        else
            fprintf('Could not find cluster profile: %s\n',mls_fname);
            cluster_profile = {};
        end
    end
end
    
%% Set up cluster properties
c = parcluster;
jobdir = fullfile(c.JobStorageLocation,sprintf('pipeline_local_%u',char(datetime('now','Format','yyyyMMddHHmmss'))));
mkdir(jobdir);
c.JobStorageLocation = jobdir;

%% Loop over cases for local processing
ncases = numel(cases);
procdir = cell(ncases,1);
for i = 1:ncases
    
    basename = sprintf('%s_%s',cases(i).UMlabel,cases(i).StudyDate);
    procdir{i} = fullfile(opts.sv_path,basename);
    fn_log = 
    
    if flag
        % Start pipeline in command window
        fprintf('%s - starting pipeline_local case #%u\n',basename,i);
        pipeline_local(basename,cases(i).Scans.Directory,procdir{i});
    else
        % Start pipeline batch job
        fprintf('%s - starting Main_Pipeline_sub as batch job: #',basename);
        j(i) = batch(@pipeline_local,1,[{basename},{cases(i).Scans.Directory},procdir(i)],cluster_profile{:});
        fprintf('%u\n',j(i).ID);
    end
end

% Wait for all jobs to complete
errflag = true(ncases,1);
dt = zeros(ncases,1);
for i = 1:ncases
    wait(job(i));
    if ~isempty(job(i).Tasks(1).Error) || strcmp(job(i).State,'failed')
        errflag(i) = false;
        getReport(job(i).Tasks(1).Error)
    else
        % Compile statistics:
        dt(i) = minutes(job(i).FinishDateTime - job(i).StartDateTime);
        val = job(i).fetchOutputs;
        val = val{1};
        if size(val,1)>1
            val = lobeTable2struct(val);
        end
        if isstruct(val)
            val = struct2table(val);
        end
        T = [T;val];
    end
    fprintf('Job %u finished after %.1f minutes.\n',i,dt(i));
    fprintf('  State: %s\n',job(i).State);
    job(i).diary

    % Delete job files:
    job(i).delete;

end
rmdir(jobdir,'s');

%% Generate files and copy command for Great Lakes execution



function GLinputs = pipeline_local(basename,expfname,insfname,procdir)

varnames = {'ID','Exp_DICOM','Ins_DICOM','ProcDir','ElxDir','ROI'};
Nv = numel(varnames);
regionnames = {'WholeLung','RL','LL','RUL','RML','RULplus','RLL','LUL','LULplus','LLi','LLL'};
Nr = numel(regionnames);
res = table('Size',[Nr,Nv],'VariableTypes',repmat({'cellstr'},1,Nv),'VariableNames',varnames,'RowNames',regionnames);
res.ROI = regionnames';

if nargin<5
    quickreg = false;
end

tt = tic;

%% Initialize folders and log file
if ~isfolder(procdir)
    mkdir(procdir);
end
fn_log = fullfile(procdir,'pipeline_log.txt');

%% Initialize parameters and results struct
fn_ext = '.nii.gz';
check_EI = false;
img = struct('flag',[false,false],'mat',{[],[]},'info',{[],[]},'label',{[],[]});
% tok = regexp(expfname,'\\([^\\\.]+)\.','tokens');
% res.ID = extractBefore(tok{1}{1},'_');
res.ProcDir(:) = {procdir};

%% Load CT images
writeLog(fn_log,'Pipeline Processing: %s\n',basename);
% Exp
writeLog(fn_log,'Loading EXP image...\n');
if isfolder(expfname)
    [img(1).mat,label,fov,orient,info] = readDICOM(expfname,[],true);
elseif ~isempty(expfname)
    [img(1).mat,label,fov,orient,info] = cmi_load(1,[],expfname);
else
    writeLog(fn_log,'  EXP image not loaded.\n');
end
img(1).flag = ~isempty(img(1).mat);
if img(1).flag
    d = size(img(1).mat);
    img(1).info = struct('label',label,'fov',fov,'orient',orient,'d',d,'voxsz',fov./d,'natinfo',info);
    img(1).info.voxvol = prod(img(1).info.voxsz);
    img(1).info.name = img(1).info.label;
end

% Ins
writeLog(fn_log,'Loading INSP image...\n');
if isfolder(insfname)
    [img(2).mat,label,fov,orient,info] = readDICOM(insfname,[],true);
elseif ~isempty(insfname)
    [img(2).mat,label,fov,orient,info] = cmi_load(1,[],insfname);
else
    writeLog(fn_log,'  INS image not loaded.\n');
end
img(2).flag = ~isempty(img(2).mat);
if img(2).flag
    d = size(img(2).mat);
    img(2).info = struct('label',label,'fov',fov,'orient',orient,'d',d,'voxsz',fov./d,'natinfo',info);
    img(2).info.voxvol = prod(img(2).info.voxsz);
    img(2).info.name = img(2).info.label;
end

res.ID(:) = {basename};
if isfolder(expfname)
    res.Exp_DICOM(:) = {expfname};
%     res.ID = sprintf('%s_%s',info.meta.PatientID,info.meta.StudyDate);
    check_EI = true;
% else
%     tok = regexp(expfname,'\\([^\\\.]+)\.','tokens');
%     res.ID = tok{1}{1};
%     if ismember('_',res.ID)
%         res.ID = extractBefore(tok{1}{1},'_');
%     end
end
res.ElxDir(:) = {fullfile(procdir,sprintf('elastix_%s',res.ID{1}))};
if isfolder(insfname)
    res.Ins_DICOM(:) = {insfname};
    check_EI = true;
end

%% Fix orientation of image:
svchk = false;
check_orient = 1; % just a way to skip check
tag = {'Exp','Ins'};
if check_orient == 1
    for ii = 1:2
        if img(ii).flag
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
    end
    clear BW
end

%% Identify Exp and Ins using lung volume; used for determining file name
if check_EI && img(1).flag && img(2).flag
    %   ** Need to double-check in case of mislabel
    
    %% Quick segmentation for Exp/Ins Identification:
    img(1).label = getRespiratoryOrgans(medfilt2_3(img(1).mat));
    img(2).label = getRespiratoryOrgans(medfilt2_3(img(2).mat));
    
    flag = (nnz(img(1).label)*img(1).info.voxvol) > (nnz(img(2).label)*img(2).info.voxvol) ...
        && (mean(img(1).mat(logical(img(1).label))) < mean(img(2).mat(logical(img(2).label))));
    if flag
        writeLog(fn_log,'Swapping INS/EXP due to lung volume\n');
        img = img([2,1]);
        tstr = res.Exp_DICOM;
        res.Exp_DICOM = res.Ins_DICOM;
        res.Ins_DICOM = tstr;
    end
    svchk = true;
elseif xor(img(1).flag,img(2).flag)
    svchk = true;
end
img(1).info.label = [res.ID{1},'_Exp'];
img(2).info.label = [res.ID{1},'_Ins'];
img(1).info.name = [res.ID{1},'_Exp'];
img(2).info.name = [res.ID{1},'_Ins'];

%% Save images:
if svchk
    %% Save nii.gz files using ID and Tag
    if img(1).flag
        fn_exp = fullfile(procdir,sprintf('%s',res.ID{1},'.exp',fn_ext));
        saveNIFTI(fn_exp,img(1).mat,img(1).info.label,img(1).info.fov,img(1).info.orient);
    end
    if img(2).flag
        fn_ins = fullfile(procdir,sprintf('%s',res.ID{1},'.ins',fn_ext));
        saveNIFTI(fn_ins,img(2).mat,img(2).info.label,img(2).info.fov,img(2).info.orient);
    end
end

%% Generate Lung Segmentation [This is when VOI don't exist]
fn_label = cell(2,1);
for itag = 1:2
    if img(itag).flag
        fn_label{itag} = fullfile(procdir,sprintf('%s.%s.label%s',res.ID{1},lower(tag{itag}),fn_ext));
        writeLog(fn_log,'%s segmentation ... ',tag{itag});
        if exist(fn_label{itag},'file')
            writeLog(fn_log,'from file\n');
            img(itag).label = readNIFTI(fn_label{itag});
        else
            img(itag).label = CTlung_Segmentation(4,img(itag).mat,img(itag).info,img(itag).info.label,procdir,fn_log);
            
%             writeLog(fn_log,'generating using YACTA\n');
%             ydir = fullfile(procdir,['yacta_',img(itag).info.name]);
%             if ~isfolder(ydir)
%                 mkdir(ydir);
%             end
%             tname = fullfile(ydir,sprintf('%s.mhd',img(itag).info.label));
%             saveMHD(tname,img(itag).mat,img(itag).info.label,img(itag).info.fov,img(itag).info.orient);
%             yacta(tname,'wait');
%             tname = dir(sprintf('%s*lung_lobes*explabels.mhd',tname));
%             img(itag).label = cmi_load(1,[],fullfile(ydir,tname(end).name));

            saveNIFTI(fn_label{itag},img(itag).label,img(itag).info.label,img(itag).info.fov,img(itag).info.orient);
        end
    end
end

%% Airways
for itag = 1:2
    if img(itag).flag
        ydir = fullfile(procdir,['yacta_',img(itag).info.name]);
        airway_res = readYACTAairways(ydir);
        if ~isempty(airway_res)
            
            res = addTableVarVal(res,'WallPct_3_8','WholeLung',airway_res.Wall_pct__3_8_);
            
            res = addTableVarVal(res,'Wall_pct',{'WholeLung','RUL','RML','RULplus','RLL','LUL','LULplus','LLi','LLL'},...
                [ airway_res.Wall_pct;...
                  airway_res.Wall_pct_RightUpperLobe;...
                  airway_res.Wall_pct_RightMidLobe;...
                  airway_res.Wall_pct_RightUpperLobePlus;...
                  airway_res.Wall_pct_RightLowerLobe;...
                  airway_res.Wall_pct_LeftUpperLobe;...
                  airway_res.Wall_pct_LeftUpperLobePlus;...
                  airway_res.Wall_pct_LeftLingula;...
                  airway_res.Wall_pct_LeftLowerLobe ]);

            res = addTableVarVal(res,'Pi10',{'WholeLung','RUL','RML','RULplus','RLL','LUL','LULplus','LLi','LLL'},...
                [ airway_res.Pi10;...
                  airway_res.Pi10_RUL;...
                  airway_res.Pi10_RML;...
                  airway_res.Pi10_RULplus;...
                  airway_res.Pi10_RLL;...
                  airway_res.Pi10_LUL;...
                  airway_res.Pi10_LULplus;...
                  airway_res.Pi10_LLi;...
                  airway_res.Pi10_LLL ]);

            res = addTableVarVal(res,'Pi15','WholeLung',airway_res.Pi15);

            genstr = {'WT'};
            segstr = { 'seg',   4   ;...
                       'subseg',5:7 };
            lobestr =  {'Right','Left','RUL','RML','RUL+','RLL','LUL','LUL+','LLi','LLL'};
            lobestr2 = {'RL','LL','RUL','RML','RULplus','RLL','LUL','LULplus','LLi','LLL'};
            for i = 1:numel(genstr)
                for j = 1:2
                    for k = 1:numel(lobestr)
                        segind = segstr{j,2};
                        vals = airway_res.(genstr{i}).(lobestr{k});
                        vals = vals(segind(segind<=numel(vals)));
                        res = addTableVarVal(res,[genstr{i},'_',segstr{j,1}],lobestr2{k},mean(vals(vals>0)));
                    end
                end
            end

            res = addTableVarVal(res,'BEI',{'WholeLung','RL','LL','RUL','RML','RULplus','RLL','LUL','LULplus','LLi','LLL'},...
                [ airway_res.BEI_Lung;...
                  airway_res.BEI_Right;...
                  airway_res.BEI_Left;...
                  airway_res.BEI_RUL;...
                  airway_res.BEI_RML;...
                  airway_res.BEI_RULplus;...
                  airway_res.BEI_RLL;...
                  airway_res.BEI_LUL;...
                  airway_res.BEI_LULplus;...
                  airway_res.BEI_LLi;...
                  airway_res.BEI_LLL ]);
        end
    end
end

%% QC segmentation
if img(1).flag
    writeLog(fn_log,'Generating EXP montage...\n');
    ind = 10:10:img(1).info.d(3);
    QCmontage('seg',cat(4,img(1).mat(:,:,ind),logical(img(1).label(:,:,ind))),img(1).info.voxsz,...
        fullfile(procdir,sprintf('%s_Montage',img(1).info.label)));
end
if img(2).flag
    writeLog(fn_log,'Generating INSP montage...\n');
    ind = 10:10:img(2).info.d(3);
    QCmontage('seg',cat(4,img(2).mat(:,:,ind),logical(img(2).label(:,:,ind))),img(2).info.voxsz,...
        fullfile(procdir,sprintf('%s_Montage',img(2).info.label)));
end

%% ScatterNet for AT on Exp CT scan
if img(1).flag
    fn_scatnet = fullfile(procdir,sprintf('%s.%s%s',res.ID{1},'scatnet',fn_ext));
    writeLog(fn_log,'Air trapping map ... ');
    if exist(fn_scatnet,'file')
        writeLog(fn_log,'from file\n');
        atMap = cmi_load(1,img(1).info.d(1:3),fn_scatnet);
    else
        writeLog(fn_log,'generating with ScatNet\n');
        atMap = ScatNet(img(1).mat,logical(img(1).label),0);
        cmi_save(0,atMap,'ScatNet',img(1).info.fov,img(1).info.orient,fn_scatnet);
    end
end

%% Vessel analysis
skip_Vessel = 1;
if skip_Vessel == 0
    if img(2).flag
        writeLog(fn_log,'Vessel analysis ...\n');
        fn_re_ins = fullfile(procdir,sprintf('re_%s.ct.nii.gz',res.ID{1}));
        fn_re_seg = fullfile(procdir,sprintf('re_%s.lobe_segmentation.nii.gz',res.ID{1}));
        if exist(fn_re_ins,'file') && exist(fn_re_seg,'file')
            tinfo = niftiinfo(fn_re_ins);
            tins = niftiread(tinfo);
            tseg = niftiread(fn_re_seg);
        else
            tinfo = init_niftiinfo(res.ID{1},img(2).info.voxsz,class(img(2).mat),img(2).info.d);
            tins = img(2).mat;
            tseg = img(2).label;
        end
        T = vesselSeg_BH( tins , tseg , tinfo , procdir );
        res = addTableVarVal(res,T);
    end
end

%% Quantify unregistered CT scans
writeLog(fn_log,'Quantifying unregistered statistics\n');
if img(1).flag
    T = CTlung_Unreg('exp',img(1).mat,img(1).info.voxvol,img(1).label,atMap);
    clear atMap;
    res = addTableVarVal(res,T);
end
if img(2).flag
    T = CTlung_Unreg('ins',img(2).mat,img(2).info.voxvol,img(2).label);
    res = addTableVarVal(res,T);
end

%% Register I2E
if img(1).flag && img(2).flag
    fn_reg = fullfile(procdir,sprintf('%s.%s%s',res.ID{1},'ins.reg',fn_ext));
    if exist(fn_reg,'file')
        writeLog(fn_log,'Loading registered INS from file...\n');
        ins_reg = readNIFTI(fn_reg);
    else
        writeLog(fn_log,'Performing registration...\n')
        lungreg_BH( img(1).mat, img(1).info, logical(img(1).label),...
                    img(2).mat, img(2).info, logical(img(2).label),...
                    res.ElxDir{1}, res.ID{1}, false, quickreg);
        % Move registered file out to procdir
        [~,outfn] = fileparts(img(2).info.name);
        outfn = fullfile(res.ElxDir{1},sprintf('%s_R.nii',outfn));
        % Elastix won't save .nii.gz, so need to resave
        ins_reg = readNIFTI(outfn);
        cmi_save(0,ins_reg,'Ins_R',img(1).info.fov,img(1).info.orient,fn_reg);
    end

    %% QC registration
    writeLog(fn_log,'Saving Registration Montage ...\n');
    ind = 10:10:img(1).info.d(3);
    QCmontage('reg',cat(4,ins_reg(:,:,ind),img(1).label(:,:,ind)),img(1).info.voxsz,...
        fullfile(procdir,sprintf('%s_Reg_Montage',res.ID{1})));
    img(2) = [];


    %% PRM calculation
    fn_PRM = fullfile(procdir,sprintf('%s.%s%s',res.ID{1},'prm',fn_ext));
    if exist(fn_PRM,'file')
        writeLog(fn_log,'Loading PRM from file...\n');
        prm10 = int8(readNIFTI(fn_PRM));
    else
        writeLog(fn_log,'Calculating PRM...\n');
        [prm10,~] = pipeline_PRM(img(1).mat,img(1).info,logical(img(1).label),ins_reg,...
            fullfile(procdir,sprintf('%s_PRM_Scatter',res.ID{1})));

        % Save PRM
        writeLog(fn_log,'Saving PRM as NIFTI ... ');
        stat = cmi_save(0,prm10,{'PRM'},img.info.fov,img.info.orient,fn_PRM);
        if stat
            writeLog(fn_log,'  PRM saved\n');
        else
            writeLog(fn_log,'  Could not save PRM to file.\n');
        end
    end
    clear ins_reg;
    
    %% Tabulate 10-color PRM results
    writeLog(fn_log,'Tabulating 10-color PRM results...\n');
    T = lobeLoop(img.label,@(mask,prm,flag)tabulatePRM(mask,prm,flag),prm10,1);
    res = addTableVarVal(res,T);

    %% map full PRM values (1:10) to (norm,fsad,emph,pd,ns)
    prmlabel = {'Norm', 'fSAD', 'Emph', 'PD',       'NS';...
                [1,2],  3,      [4,5],  [8,9,10],   6    };
    prm5 = prm10;
    for i = 1:size(prmlabel,2)
        prm5(ismember(prm10,prmlabel{2,i})) = i;
    end
    clear prm10

    %% QC PRM
    writeLog(fn_log,'Generating PRM Montage ...\n');
    QCmontage('prm',cat(4,img.mat(:,:,ind),double(prm5(:,:,ind))),...
        img.info.voxsz,fullfile(procdir,sprintf('%s_PRM_Montage',res.ID{1})));

    %% Tabulate 5-color PRM results
    writeLog(fn_log,'Tabulating 5-color PRM results...\n');
    T = lobeLoop(img.label,@(mask,prm,flag)tabulatePRM(mask,prm,flag),prm5,0);
    res = addTableVarVal(res,T);

    %% Calculate tPRM
    prmlabel = ["norm","fsad","emph","pd"];
    mflabel = ["v","s","b","x"];
    fn_tprm = fullfile(procdir,...
        string(res.ID{1})+".tprm."+prmlabel+"."+mflabel'+string(fn_ext));
    if all(cellfun(@(x)exist(x,'file'),fn_tprm))
        writeLog(fn_log,'Loading tPRM from files ...\n');
        clear prm5 prm10;
        for iprm = 1:numel(prmlabel)
            for imf = 1:numel(mflabel)
                writeLog(fn_log,'   %s - %s\n',prmlabel(iprm),mflabel(imf))
                tprm = readNIFTI(fn_tprm(imf,iprm));
                str = prmlabel(iprm)+'_'+upper(mflabel(imf));
                T = lobeLoop(img.label,@(mask,tprm,str)tabulateTPRM(mask,tprm,str),tprm,str);
                res = addTableVarVal(res,T);
            end
        end
    else
        t = tic;
        writeLog(fn_log,'Generating tPRM ...\n');

        %% Calculate MF values
        p = minkowskiFun(prm5,'thresh',1:4,...
                             'tmode','==',...
                             'n',10*ones(1,3),...
                             'gridsp',5,...
                             'voxsz',img.info.voxsz,...
                             'mask',logical(img.label),...
                             'prog',0);
        clear prm5;

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
                str = prmlabel(ithresh)+'_'+upper(mflabel(imf));
                T = lobeLoop(img.label,@(mask,tprm,str)tabulateTPRM(mask,tprm,str),tprm,str);
                res = addTableVarVal(res,T);
            end
        end
        writeLog(fn_log,'... tPRM complete (%s)\n',datestr(duration(0,0,toc(t)),'HH:MM:SS'));

    end
end
    
%% Save Results Table:
writetable(res,fullfile(procdir,[res.ID{1},'_PipelineResults.csv']));

writeLog(fn_log,'Pipeline total time = %s\n',datestr(duration(0,0,toc(tt)),'HH:MM:SS'))


function T = tabulatePRM(mask,prm,flag)
    if flag % 10-color
        vals = 1:10;
        tag = cellfun(@num2str,num2cell(vals),'UniformOutput',false);
    else % 5-color
        vals = 1:5;
        tag = {'Norm', 'fSAD', 'Emph', 'PD', 'NS'};
    end
    nv = numel(vals);
    vname = strcat('PRM_',tag);
    T = table('Size',[1,nv],'VariableTypes',repmat({'double'},1,nv),...
        'VariableNames',vname);
    np = nnz(mask);
    for i = 1:numel(vals)
        T.(vname{i}) = nnz(prm(mask)==i)/np*100;
    end
    
function T = tabulateTPRM(mask,tprm,str)
    vname = sprintf('tPRM_%s',str);
    T = table('Size',[1,1],'VariableTypes',{'double'},'VariableNames',{vname});
    T.(vname) = mean(tprm(mask));
    
function img = medfilt2_3(img)
    for i = 1:size(img,3)
        img(:,:,i) = medfilt2(img(:,:,i));
    end
    
function T = addTableVarVal(T,varargin)
    if nargin==2 && istable(varargin{1}) && strcmp(varargin{1}.Properties.VariableNames{1},'LOBE')
        vals = varargin{1};
        regionstr = vals.LOBE;
        vals = removevars(vals,'LOBE');
        varstr = vals.Properties.VariableNames;
    elseif nargin==4
        varstr = varargin{1};
        regionstr = varargin{2};
        vals = varargin{3};
        if ischar(regionstr)
            regionstr = {regionstr};
        end
        if ischar(varstr)
            varstr = {varstr};
        end
        if numel(varstr)<size(vals,2)
            varstr = strcat(varstr,'_',cellfun(@num2str,num2cell(1:size(vals,2)),'UniformOutput',false));
        end
    else
        return;
    end
    if numel(regionstr)==size(vals,1)
        for i = 1:numel(varstr)
            if istable(vals)
                tvals = vals.(vals.Properties.VariableNames{i});
            else
                tvals = vals(:,i);
            end
            switch class(tvals)
                case 'cell'
                    defval = {''};
                case 'double'
                    defval = nan;
                otherwise
                    defval = {};
            end
            if ~ismember(varstr{i},T.Properties.VariableNames)
                T = addvars(T,repmat(defval,size(T,1),1),'NewVariableNames',varstr(i));
            end
            T.(varstr{i})(regionstr) = tvals;
        end
    end
    