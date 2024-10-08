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
opts = struct('username','',...  % UMID
    'dcm_path','',...  % Location of DICOM data
    'save_path','',... % Directory containing case folders (procdirs)
    'quickreg',false);

% if nargin
%     batchtype = find(strcmp(batchtype,{'serial','batch','serialbatch'}),1);
%     if isempty(batchtype)
%         error('Invalid batchtype input.');
%     end
% else % Default to 'serial'
%     batchtype = 1;
% end

% Find files
[cases,opts] = catalog_select_3('opts',opts);

% Make sure the save directory is on Turbo, for access from GL
[~,flag_turbo] = checkTurboPath(opts.save_path);
if flag_turbo
    error('Save path must be on Turbo for access from Great Lakes.');
end

% % Determine batch setup
% if batchtype==3
%     cluster_profile = {'Profile','cmi_1'};
%     if ~ismember('cmi_1',parallel.clusterProfiles)
%         mls_fname = fullfile(fileparts(which('cmi')),'cmi_1.mlsettings');
%         if exist(mls_fname,'file')
%             fprintf('Importing cluster profile: %s\n',mls_fname);
%             parallel.importProfile(mls_fname);
%         else
%             fprintf('Could not find cluster profile: %s\n',mls_fname);
%             cluster_profile = {};
%         end
%     end
% end

% Set up cluster properties
if batchtype~=1
    c = parcluster;
    jobdir = fullfile(c.JobStorageLocation,sprintf('pipeline_local_%u',char(datetime('now','Format','yyyyMMddHHmmss'))));
    mkdir(jobdir);
    c.JobStorageLocation = jobdir;
end

% Run as batch with parpool size of 6
batch(@pipeline_loop,0,{},'Pool',6);


% Loop over cases for local processing
% ncases = numel(cases);
% procdir = cell(ncases,1);
for i = 1:ncases
    basename = sprintf('%s_%s',cases(i).UMlabel,cases(i).StudyDate);
    procdir{i} = fullfile(opts.save_path,basename);
    if batchtype==1
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
if batchtype~=1
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
if batchtype~=1
    rmdir(jobdir,'s');
end

GL_run(opts.username, 'Main_Pipeline_GL_sub', {procdir,opts}, [true,false], [false,true],...
    'ProcessMemory',24,'ProcessTime',60*12)

function procdir = pipeline_loop(cases,save_path)
% Loop function for inside batch process
ncases = numel(cases);
procdir = cell(ncases,1);
parfor i = 1:ncases
    basename = sprintf('%s_%s',cases(i).UMlabel,cases(i).StudyDate);
    procdir{i} = fullfile(save_path,basename);
    pipeline_local(basename,cases(i).Scans.Directory,procdir{i});
end

function res = pipeline_local(ID,expfname,insfname,procdir)
% Local execution for YACTA segmentation and airways analysis
% * procdir must be on Turbo for Great Lakes to have access to it

t = tic;

% Initialize results table
varnames = {'ID','Exp_DICOM','Ins_DICOM','ProcDir'};
Nv = numel(varnames);
res = table('Size',[1,Nv],'VariableTypes',repmat({'cellstr'},1,Nv),'VariableNames',varnames);
res.ID = {ID};
res.ProcDir = procdir;

% Initialize folders and log file
if ~isfolder(procdir)
    mkdir(procdir);
end
fn_log = fullfile(procdir,'pipeline_log.txt');
writeLog(fn_log,'%s\n',ID);

% Initialize parameters and results struct
fn_ext = '.nii.gz';
img = struct('flag',[false,false],'mat',{[],[]},'info',{[],[]},'label',{'',''},'dcmpath',{'',''});

fn = fullfile(procdir,string(ID)+[".exp",".exp.label";".ins",".ins.label"]+fn_ext);
fnflag = cellfun(@(x)exist(x,'file'),fn);
if any(~fnflag(:,1))
    % If either image is missing, must redo all from DICOM
    fnflag(:) = false;
end

% Load CT data
tagstr = {'Exp','Ins'};
fn_in = {expfname,insfname};
for i = 1:2
    orientchk = false;
    % CT from DICOM
    writeLog(fn_log,'%s : CT ... ',tagstr{i});
    if fnflag(i,1)
        writeLog(fn_log,'file found: %s\n',fn{i,1});
        if ~fnflag(i,2)
            % Need to load image for segmentation
            writeLog(fn_log,'   Loading for segmentation.');
            [img(i).mat,label,fov,orient,info] = cmi_load(1,[],fn_in{i});
        end
    else
        writeLog(fn_log,'   loading data ... ',fn{i,1});
        if isfolder(fn_in{i})
            writeLog(fn_log,'from DICOM');
            [img(i).mat,label,fov,orient,info] = readDICOM(fn_in{i},[],true);
            img(i).dcmpath = fn_in{i};
            orientchk = true;
        elseif ~isempty(fn_in{i})
            writeLog(fn_log,'from file')
            [img(i).mat,label,fov,orient,info] = cmi_load(1,[],fn_in{i});
        else
            writeLog(fn_log,'image not loaded.\n');
        end
    end
    
    % Set image info
    img(i).flag = ~isempty(img(i).mat);
    if img(i).flag
        d = size(img(i).mat);
        img(i).info = struct('label',label,'fov',fov,'orient',orient,'d',d,'voxsz',fov./d,'natinfo',info);
        img(i).info.voxvol = prod(img(i).info.voxsz);
        img(i).info.name = img(i).info.label;
    end
    
    % Check orientation using a bone threshold
    if orientchk
        % Find orientation of shoulder bones to see if permute is needed
        BW = img(i).mat(:,:,end) > -150;
        %BW = max(img(ii).mat(:,:,round(img(ii).info.d(3)/2):end)>800,[],3);
        prop = regionprops(BW,'Orientation','Area');
        if mod(round(prop([prop.Area]==max([prop.Area])).Orientation/90),2)
            writeLog(fn_log,'Permuting %s\n',tagstr{i});
            img(i).mat = permute(img(i).mat,[2,1,3]);
            img(i).info.voxsz = img(i).info.voxsz([2,1,3]);
            img(i).info.fov = img(i).info.fov([2,1,3]);
            img(i).info.d = img(i).info.d([2,1,3]);
            img(i).info.orient = img(i).info.orient([2,1,3,4],[2,1,3,4]);
        end
    end
end

% Check for Exp/Ins
if any(~fnflag(:,1)) && img(1).flag && img(2).flag
    %   ** Need to double-check in case of mislabel
    
    % Quick segmentation for Exp/Ins Identification:
    seg1 = getRespiratoryOrgans(medfilt2_3(img(1).mat));
    seg2 = getRespiratoryOrgans(medfilt2_3(img(2).mat));
    
    flag = (nnz(seg1)*img(1).info.voxvol) > (nnz(seg2)*img(2).info.voxvol) ...
        && (mean(img(1).mat(logical(seg1))) < mean(img(2).mat(logical(seg2))));
    if flag
        writeLog(fn_log,'Swapping INS/EXP due to lung volume\n');
        img = img([2,1]);
    end
    clear seg1 seg2
end
res.Exp_DICOM = {img(1).dcmpath};
res.Ins_DICOM = {img(2).dcmpath};
img(1).info.label = [ID,'_Exp'];
img(2).info.label = [ID,'_Ins'];
img(1).info.name =  [ID,'_Exp'];
img(2).info.name =  [ID,'_Ins'];

% Generate YACTA segmentations
for i = 1:2
    writeLog(fn_log,'%s : Segmentation ... ',tagstr{i});
    if fnflag(i,2)
        writeLog(fn_log,'file found: %s\n',fn{i,2});
    else
        writeLog(fn_log,'generating new ...\n');
        seg = CTlung_Segmentation(4,img(i).mat,img(i).info,img(i).info.label,procdir,fn_log);
    end
    
    if img(i).flag
        % Save image:
        if ~fnflag(i,1)
            saveNIFTI(fn{i,1},img(i).mat,img(i).info.label,img(i).info.fov,img(i).info.orient);
        end
        % Save segmentation
        if ~fnflag(i,2)
            saveNIFTI(fn{i,2},seg,img(i).info.label,img(i).info.fov,img(i).info.orient);
        end
    end
end

dt = toc(t);
writeLog(fn_log,'Local processing completed after %.1f minutes\n',dt/60);

% Save results to CSV file:
writetable(res,fullfile(procdir,[ID,'_PipelineResults.csv']));

function img = medfilt2_3(img)
for i = 1:size(img,3)
    img(:,:,i) = medfilt2(img(:,:,i));
end



