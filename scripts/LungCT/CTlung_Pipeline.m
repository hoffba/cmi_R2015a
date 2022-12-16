function results = CTlung_Pipeline(varargin)

results = [];

try
    
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
[cases,opts] = catalog_select_4('opts',opts);
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
    cases(i).fname_exp = cases(i).Scans(1).DataPath;
    cases(i).fname_ins = cases(i).Scans(2).DataPath;
    cases(i).basename = sprintf('%s_%s',cases(i).UMlabel,cases(i).StudyDate);
    cases(i).procdir = fullfile(opts.save_path,cases(i).basename);
end

% Determine how to run the pipeline, based on opts.cluster:
%  - 'GL'    = run DICOM loading and YACTA locally, then the rest on GL
%  - 'batch' = run full pipeline locally in batch job
%  - 'debug' = run full pipeline locally in command window
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
    c = parcluster('local');
    nworkers_orig = c.NumWorkers;
    jobdir = fullfile(c.JobStorageLocation,sprintf('pipeline_local_%s',opts.timestamp));
    mkdir(jobdir);
    c.JobStorageLocation = jobdir;
    np = min([c.NumWorkers-1,ncases,opts.par_size]);
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

    
catch err
    disp('check');
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
                cases(i).basename,cases(i).Scans.DataPath,cases(i).procdir,opts);
        end
        % Flag processes as the complete
        for i = 1:ncases
            idx = fetchNext(f);
            fprintf('Finished process #%d of %d: %s\n', idx, ncases, cases(idx).basename);
        end
    otherwise % 'batch' or 'debug' - Run full process locally
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
CTlung_Pipeline_local(case_i.basename,case_i.Scans.DataPath,case_i.procdir,opts);
res = CTlung_Pipeline_sub(case_i.procdir,opts);