function Main_Pipeline_GL(batchtype)
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
opts = struct('dcmpath','',...
              'sv_path','',...
              'quickreg',false);

if nargin
    batchtype = find(strcmp(batchtype,{'serial','batch','serialbatch'}),1);
    if isempty(batchtype)
        error('Invalid batchtype input.');
    end
else % Default to 'serial'
    batchtype = 1;
end
    
% Find files
[cases,opts] = catalog_select_3(opts);

% Make sure the save directory is on Turbo, for access from GL
[~,flag_turbo] = checkTurboPath(opts.savepath);
if flag_turbo
    error('Save path must be on Turbo for access from Great Lakes.');
end

% Determine batch setup
if batchtype==3
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

% Set up cluster properties
if ~flag
    c = parcluster;
    jobdir = fullfile(c.JobStorageLocation,sprintf('pipeline_local_%u',char(datetime('now','Format','yyyyMMddHHmmss'))));
    mkdir(jobdir);
    c.JobStorageLocation = jobdir;
end

% Loop over cases for local processing
ncases = numel(cases);
procdir = cell(ncases,1);
for i = 1:ncases
    basename = sprintf('%s_%s',cases(i).UMlabel,cases(i).StudyDate);
    procdir{i} = fullfile(opts.sv_path,basename);
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
if ~flag
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
        end
        fprintf('Job %u finished after %.1f minutes.\n',i,dt(i));
        fprintf('  State: %s\n',job(i).State);
        job(i).diary

        % Delete job files:
        job(i).delete;

    end
end
rmdir(jobdir,'s');

GL_run(username, 'Main_Pipeline_GL_sub', {procdir,opts}, [true,false], [false,true])


function res = pipeline_local(ID,expfname,insfname,procdir)
% Local execution for YACTA segmentation and airways analysis
% * procdir must be on Turbo for Great Lakes to have access to it

    t = tic;

    % Initialize results table
    varnames = {'ID','Exp_DICOM','Ins_DICOM','ProcDir'};
    Nv = numel(varnames);
    res = table('Size',[1,Nv],'VariableTypes',repmat({'cellstr'},1,Nv),'VariableNames',varnames);
    res.ID = {ID};
    
    % Initialize folders and log file
    if ~isfolder(procdir)
        mkdir(procdir);
    end
    fn_log = fullfile(procdir,'pipeline_log.txt');

    % Initialize parameters and results struct
    fn_ext = '.nii.gz';
    check_EI = false;
    img = struct('flag',[false,false],'mat',{[],[]},'info',{[],[]},'label',{[],[]});

    % Load CT images
    writeLog(fn_log,'Pipeline Processing: %s\n',ID);
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

    if isfolder(expfname)
        res.Exp_DICOM = {expfname};
        check_EI = true;
    end
    if isfolder(insfname)
        res.Ins_DICOM = {insfname};
        check_EI = true;
    end

    % Fix orientation of image:
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

    % Identify Exp and Ins using lung volume; used for determining file name
    if check_EI && img(1).flag && img(2).flag
        %   ** Need to double-check in case of mislabel

        % Quick segmentation for Exp/Ins Identification:
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
    img(1).info.label = [ID,'_Exp'];
    img(2).info.label = [ID,'_Ins'];
    img(1).info.name =  [ID,'_Exp'];
    img(2).info.name =  [ID,'_Ins'];

    % Save images:
    if svchk
        % Save nii.gz files using ID and Tag
        if img(1).flag
            fn_exp = fullfile(procdir,sprintf('%s',ID,'.exp',fn_ext));
            saveNIFTI(fn_exp,img(1).mat,img(1).info.label,img(1).info.fov,img(1).info.orient);
        end
        if img(2).flag
            fn_ins = fullfile(procdir,sprintf('%s',ID,'.ins',fn_ext));
            saveNIFTI(fn_ins,img(2).mat,img(2).info.label,img(2).info.fov,img(2).info.orient);
        end
    end

    % Generate Lung Segmentation [This is when VOI don't exist]
    fn_label = cell(2,1);
    for itag = 1:2
        if img(itag).flag
            fn_label{itag} = fullfile(procdir,sprintf('%s.%s.label%s',ID,lower(tag{itag}),fn_ext));
            writeLog(fn_log,'%s segmentation ... ',tag{itag});
            if exist(fn_label{itag},'file')
                writeLog(fn_log,'from file\n');
                img(itag).label = readNIFTI(fn_label{itag});
            else
                img(itag).label = CTlung_Segmentation(4,img(itag).mat,img(itag).info,img(itag).info.label,procdir,fn_log);
                saveNIFTI(fn_label{itag},img(itag).label,img(itag).info.label,img(itag).info.fov,img(itag).info.orient);
            end
        end
    end

    dt = toc(t);
    writeLog(fn_log,'Local processing completed after %.1f minutes\n',dt/60);
    
    % Save results to CSV file:
    writetable(res,fullfile(procdir,[ID,'_PipelineResults.csv']));
    
    


