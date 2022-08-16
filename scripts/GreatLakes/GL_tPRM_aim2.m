function GL_tPRM_aim2(x,varargin)


if ischar(x) %%%% Local Execution %%%%
    username = x;
    GLstr = 'GL_tPRM_aim2';
    
    % Parse inputs
    N = [];
    if nargin==1
        selpath = uigetdir(pwd,'Select Home directory for 5yr Reg Data');
        [~,shortpath,flag,LOC_turbo_path] = checkTurboPath(selpath);
        if flag
            warning('Selected path is not on Turbo storage.');
            return;
        end
        % Find Subject directories in selected path
        fn = dir(selpath);
        fn(contains({fn.name},{'output'})) = [];
        fn(contains({fn.name},{'.'})) = [];
        fn = {fn.name}';
        shortpath = repmat(shortpath,numel(fn),1);
    elseif ischar(varargin{1})
        % Find Subject directories in selected path
        [~,shortpath,flag,LOC_turbo_path] = checkTurboPath(varargin{1});
        if flag
            warning('Input path is not on Turbo storage.');
            return;
        end
        fn = dir(varargin{1});
        fn(contains({fn.name},{'output'})) = [];
        fn(contains({fn.name},{'.'})) = [];
        fn = {fn.name};
        shortpath = repmat(shortpath,numel(fn),1);
    elseif iscellstr(varargin{1})
        [~,shortpath,flag,LOC_turbo_path] = checkTurboPath(varargin{1});
        if all(flag)
            warning('Inputs are not on Turbor storage.');
            return
        end
    elseif isnumeric(varargin{1})
        N = varargin{1};
    end
    nf = numel(fn);
    
    % Determine batch job inputs and save input file for GL execution:
    tstr = char(datetime('now','Format','yyyyMMddHHmmss'));
    jobname = sprintf('%s_%s',GLstr,tstr);
    fname = [jobname,'.sh'];
    tempdir = fullfile(LOC_turbo_path,'GreatLakes','temp');
    fn_inputs = fullfile(tempdir,sprintf('batch_inputs_%s.mat',tstr));
    GL_turbo_dir = {'/nfs/turbo/umms-cgalban'};
    fn = cellfun(@(x,y)strjoin([GL_turbo_dir,x(2:end),{y}],'/'),shortpath,fn,'UniformOutput',false);
    
    % Write sbatch .sh and inputs file:
    inputs = struct('fn',{fn});
    per_mem = 8; % memory estimate per case (GB)
    per_time = 15; % time estimate per case (minutes)
    GL_sbatch(fullfile(tempdir,fname),jobname,GLstr,username,fn_inputs,inputs,nf,per_mem,per_time,N);
    
    % Run SBATCH
    sb_str = sprintf('cd /nfs/turbo/umms-cgalban/GreatLakes/temp && sbatch %s',fname);
    clipboard('copy',sb_str);
    fprintf(['\n\n',...
        '1) Use PuTTy to log onto Great Lakes:\n',...
        '   Host Name: greatlakes.arc-ts.umich.edu\n',...
        '   Login: Level 1\n',...
        '2) Start SBATCH using terminal commands:\n',...
        '   (ALREADY IN CLIPBOARD, just paste into PuTTy using middle mouse button)\n',...
        '   ',sb_str,'\n'],fname);
    
else %%%% Great Lakes Execution %%%%
    
    % In GL, load file containing relevant inputs to the processing function
    %   and start processing jobs
    
    if ischar(varargin{1}) && exist(varargin{1},'file') && strcmp(varargin{1}(end-3:end),'.mat')
        inputs_fname = varargin{1};
    else
        error('Invalid input: must be valid .mat file containing function inputs.');
    end
    
    fprintf('Available workers on this node: %s\n',getenv('SLURM_CPUS_PER_TASK'));
    
    p = load(inputs_fname);
    flds = {'fn'};
    ind = ~isfield(p,flds);
    if any(ind)
        error(['Invalid input. Missing:',repmat(' %s',1,nnz(ind))],flds{ind});
    elseif ~iscellstr(p.fn)
        error('Invalid input. fn must be cellstr.');
    end
    if isfield(p,'Nnodes')
        njobs = p.Nnodes;
    else
        njobs = 1;
    end
    
    % Check for job array:
    slurm = getSLURM;
    jobstr = slurm.jobname;
    if ~isempty(slurm.arraytaskid)
        fprintf('Running job number %u of %u\n',slurm.arraytaskid,njobs);
        jobstr = sprintf('%s_array%u',jobstr,slurm.arraytaskid);
    end
    nf = numel(p.fn);
    if ~isempty(slurm.arraytaskid)
        ind = round((nf/njobs)*[slurm.arraytaskid-1 slurm.arraytaskid]);
        
        fprintf('ind: %u-%u\n',ind+[1,0]);
        
        ind = (ind(1)+1):ind(2);
    else
        ind = 1:nf;
    end
    nf = numel(ind);
    
    % Set up cluster properties
    c = parcluster;
    jobdir = fullfile(c.JobStorageLocation,jobstr);
    if ~isfolder(jobdir)
        mkdir(jobdir);
    end
    c.JobStorageLocation = jobdir;
    
    parpool(c,slurm.cpupertask-1);
    
    T = [];
    procdir = p.fn;
    parfor i = 1:numel(procdir)
        tT = tPRM_aim2_sub(procdir{i});
        T = [T;tT];
    end
    save(fullfile(fileparts(procdir{1}),sprintf('tPRM_aim2_Results_%s',jobstr)),'T');
end