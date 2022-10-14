function Main_Pipeline_GL_2(opts,debug)
% Start CT lung pipeline script (Main_Pipeline_sub) with GUI selection for
% cases based on DCMcatalog.csv

% cluster_profile = {};
if ~nargin
    opts = struct('username','',...  % UMID
        'dcm_path','',...  % Location of DICOM data
        'save_path','',... % Directory containing case folders (procdirs)
        'quickreg',false);
end
if nargin<2
    debug = false;
end

% Find files
[cases,opts] = catalog_select_3('opts',opts);
ncases = numel(cases);

% Make sure the save directory is on Turbo, for access from GL
[~,flag_turbo] = checkTurboPath(opts.save_path);
if flag_turbo
    error('Save path must be on Turbo for access from Great Lakes.');
end

for i = 1:ncases
    cases(i).basename = sprintf('%s_%s',cases(i).UMlabel,cases(i).StudyDate);
    cases(i).procdir = fullfile(opts.save_path,cases(i).basename);
end

if debug
    for i = 1:ncases
        pipeline_local(cases(i).basename,cases(i).Scans.Directory,cases(i).procdir,opts);
    end
else
    % Set up cluster properties
    c = parcluster;
    jobdir = fullfile(c.JobStorageLocation,sprintf('pipeline_local_%u',char(datetime('now','Format','yyyyMMddHHmmss'))));
    mkdir(jobdir);
    c.JobStorageLocation = jobdir;

    % Run as batch with parpool size of 6
    myCluster = parcluster('local');
    np = min(myCluster.NumWorkers-1,numel(cases));
    fprintf('Running local processes as batch with pool size of %u\nPlease wait at least an hour before starting GL process\n',np);
    batch(myCluster,@pipeline_loop,0,{cases,opts},'Pool',np);
end

GL_run(opts.username, 'Main_Pipeline_GL_sub', {{cases.procdir}',opts}, [true,false], [false,true],...
    'ProcessMemory',24,'ProcessTime',720)

function pipeline_loop(cases,opts)
% Loop function for inside batch process
ncases = numel(cases);

% Start queue for cases:
f(1:ncases) = parallel.FevalFuture;
for i = 1:ncases
    fprintf('Starting processing for case #%d of %d: %s\n',i,ncases,cases(i).basename);
    f(i) = parfeval(@(x,y,z,k,l)pipeline_local(x,y,z,k,l),0,cases(i).basename,cases(i).Scans.Directory,cases(i).procdir,opts);
end

% Flag processes as they complete
for i = 1:ncases
    completedIdx = fetchNext(f);
    fprintf('Finished process #%d of %d: %s\n', completedIdx, ncases, cases(completedIdx).basename);
end

% parfor i = 1:ncases
%     pipeline_local(cases(i).basename,cases(i).Scans.Directory,cases(i).procdir);
% end

function res = pipeline_local(ID,expfname,insfname,procdir,opts)
% Local execution for YACTA segmentation and airways analysis
% * procdir must be on Turbo for Great Lakes to have access to it

t = tic;

% Initialize results table
regionnames = {'WholeLung','RL','LL','RUL','RML','RULplus','RLL','LUL','LULplus','LLi','LLL'};
Nr = numel(regionnames);
varnames = {'ID','Exp_DICOM','Ins_DICOM'};
Nv = numel(varnames);
res = table('Size',[Nr,Nv],'VariableTypes',repmat({'cellstr'},1,Nv),'VariableNames',varnames);
res.ID(:) = {ID};

% Initialize folders and log file
if ~isfolder(procdir)
    mkdir(procdir);
end
fn_log = fullfile(procdir,sprintf('pipeline_log.txt'));
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
            writeLog(fn_log,'   Loading for segmentation.\n');
            [img(i).mat,label,fov,orient,info] = cmi_load(1,[],fn_in{i});
        end
    else
        writeLog(fn_log,'   loading data ... ',fn{i,1});
        if isfolder(fn_in{i})
            writeLog(fn_log,'from DICOM\n');
            [img(i).mat,label,fov,orient,info] = readDICOM(fn_in{i},[],true);
            img(i).dcmpath = fn_in{i};
            orientchk = true;
        elseif ~isempty(fn_in{i})
            writeLog(fn_log,'from file\n')
            [img(i).mat,label,fov,orient,info] = cmi_load(1,[],fn_in{i});
        else
            writeLog(fn_log,'image not loaded.\n');
        end
    end
    
    img(i).flag = ~isempty(img(i).mat);
    if img(i).flag
        % Set image info
        d = size(img(i).mat);
        img(i).info = struct('label',label,'fov',fov,'orient',orient,'d',d,'voxsz',fov./d,'natinfo',info);
        img(i).info.voxvol = prod(img(i).info.voxsz);
        img(i).info.name = img(i).info.label;
        
        % Check orientation using a bone threshold
        if orientchk && opts.orient_check
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
res.Exp_DICOM(:) = {img(1).dcmpath};
res.Ins_DICOM(:) = {img(2).dcmpath};
img(1).info.label = [ID,'_Exp'];
img(2).info.label = [ID,'_Ins'];
img(1).info.name =  [ID,'_Exp'];
img(2).info.name =  [ID,'_Ins'];

% Generate YACTA segmentations
for i = 1:2
    if img(i).flag
        writeLog(fn_log,'%s : Segmentation ... ',tagstr{i});
        if fnflag(i,2)
            writeLog(fn_log,'file found: %s\n',fn{i,2});
        else
            writeLog(fn_log,'generating new ...\n');
            seg = CTlung_Segmentation(4,img(i).mat,img(i).info,img(i).info.label,procdir,fn_log);
        end
        
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



