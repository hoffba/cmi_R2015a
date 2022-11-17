function results = Main_Pipeline_GL_2(opts)
% Start CT lung pipeline script (Main_Pipeline_sub) with GUI selection for
% cases based on DCMcatalog.csv

% cluster_profile = {};
if ~nargin
    opts = struct('username','',...  % UMID
        'dcm_path','',...  % Location of DICOM data
        'save_path','',... % Directory containing case folders (procdirs)
        'quickreg',false);
end

% Find files
[cases,opts] = catalog_select_3('opts',opts);
if isempty(cases)
    return;
end
ncases = numel(cases);

% Make sure the save directory is on Turbo, for access from GL
if strcmp(opts.cluster,'GL')
    [~,flag_turbo] = checkTurboPath(opts.save_path);
    if flag_turbo
        error('Save path must be on Turbo for access from Great Lakes.');
    end
end

for i = 1:ncases
    cases(i).basename = sprintf('%s_%s',cases(i).UMlabel,cases(i).StudyDate);
    cases(i).procdir = fullfile(opts.save_path,cases(i).basename);
end

% Determine how to run the pipeline, based on opts.cluster:
%  - 'GL'    = run DICOM loading and YACTA locally, then the rest on GL
%  - 'batch' = run full pipeline locally in batch job
%  - 'debug' = run full pipeline locally in command window
c = parcluster('local');
nworkers_orig = c.NumWorkers;
opts.timestamp = char(datetime('now','Format','yyyyMMddHHmmss'));
if strcmp(opts.cluster,'debug')
    
%~~~~~~~~~ debug ~~~~~~~~~
    results = [];
    for i = ncases
        res = pipeline_full(cases(i),opts);
        if ~isempty(res) && istable(res) && (isempty(results) || (size(res,2)==size(results,2)))
            results = [results;res];
        end
    end
else
    % Set up cluster properties
    jobdir = fullfile(c.JobStorageLocation,sprintf('pipeline_local_%s',opts.timestamp));
    mkdir(jobdir);
    c.JobStorageLocation = jobdir;
    np = min([c.NumWorkers-1,numel(cases),opts.par_size]);
    if strcmp(opts.cluster,'GL') % Run DICOM load and YACTA locally, then the rest on GL
        
%~~~~~~~~~ GL ~~~~~~~~~
        fprintf(['Running local processes as batch with pool size of %d\n',...
                 'Please wait at least an hour before starting GL process\n'],np);
        batch(@pipeline_loop,1,{cases,opts},'Pool',nworkers_orig,'Pool',np);
        GL_run(opts.username, 'Main_Pipeline_GL_sub', {{cases.procdir}',opts}, [true,false], [false,true],...
            'ProcessMemory',24,'ProcessTime',720,'TimeStamp',opts.timestamp)
    else
        
%~~~~~~~~~ batch ~~~~~~~~~
        fprintf('Running FULL PIPELINE in batch with pool size of %d\n',np);
        batch(@pipeline_loop,1,{cases,opts},'Pool',np);
    end
end

function results = pipeline_loop(cases,opts)
results = [];

% Loop function for inside batch process
ncases = numel(cases);

f(1:ncases) = parallel.FevalFuture;
switch opts.cluster
    case 'GL' 
        % Start queue for local processes:
        for i = 1:ncases
            fprintf('Starting processing for case #%d of %d: %s\n',i,ncases,cases(i).basename);
            f(i) = parfeval(@(x,y,z,k,l)pipeline_local(x,y,z,k,l),0,...
                cases(i).basename,cases(i).Scans.Directory,cases(i).procdir,opts);
        end
        % Flag processes as the complete
        for i = 1:ncases
            idx = fetchNext(f);
            fprintf('Finished process #%d of %d: %s\n', idx, ncases, cases(idx).basename);
        end
    case 'batch' % Run full process locally in a batch job
        % Start queue for parallel processes:
        for i = 1:ncases
            fprintf('Starting processing for case #%d of %d: %s\n',i,ncases,cases(i).basename);
            f(i) = parfeval(@(x,y)pipeline_full(x,y),1,cases(i),opts);
        end
        % Flag processes as they complete
        for i = 1:ncases
            [idx,res] = fetchNext(f);
            fprintf('Finished process #%d of %d: %s\n', idx, ncases, cases(idx).basename);
            % Gather results
            if ~isempty(res) && istable(res) &&(isempty(results) || (size(res,2)==size(results,2)))
                results = [results;res];
            end
        end
        % Write compiled results to CSV
        if ~isempty(results) && istable(results)
            writetable(results,fullfile(opts.save_path,'Pipeline_Results.csv'));
        end
end

function res = pipeline_full(case_i,opts)
pipeline_local(case_i.basename,case_i.Scans.Directory,case_i.procdir,opts);
res = Main_Pipeline_GL_sub(case_i.procdir,opts);

function res = pipeline_local(ID,expfname,insfname,procdir,opts)
% Local execution for YACTA segmentation and airways analysis
% * procdir must be on Turbo for Great Lakes to have access to it

t = tic;
try
    
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
    img = struct('flag',[false,false],'mat',{[],[]},'info',{[],[]},'label',{'',''},'dcmpath',{'Unknown','Unknown'});

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
                if size(img(i).mat,4)>1
                    img(i).mat(:,:,:,2:end) = [];
                    img(i).info(2:end) = [];
                    label(2:end) = [];
                end
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
%                 BW = img(i).mat(:,:,end) > -150;
                BW = max(img(i).mat(:,:,round(img(i).info.d(3)/2):end)>800,[],3);
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
            if ~fnflag(i,2) && numel(seg)>1
                saveNIFTI(fn{i,2},seg,img(i).info.label,img(i).info.fov,img(i).info.orient);
            end
        end
    end

    dt = toc(t);
    writeLog(fn_log,'Local processing completed after %.1f minutes\n',dt/60);

    % Save results to CSV file:
    writetable(res,fullfile(procdir,[ID,'_PipelineResults.csv']));
catch err
    writeLog(fn_log,'Pipeline ERROR:\n%s',getReport(err));
end

function img = medfilt2_3(img)
for i = 1:size(img,3)
    img(:,:,i) = medfilt2(img(:,:,i));
end



