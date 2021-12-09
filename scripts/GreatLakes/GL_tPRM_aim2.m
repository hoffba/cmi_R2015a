function GL_tPRM_aim2(varargin)

if (nargin==0) || ~isnumeric(varargin{1}) %%%% Local Execution %%%%
    
    GLstr = 'GL_tPRM_aim2';
    
    if ~nargin
        selpath = uigetdir(pwd,'Select Home directory for 5yr Reg Data');
        [fn,shortpath,flag,LOC_turbo_path] = checkTurboPath(selpath);
        if flag
            warning('Selected path is not on Turbo storage.');
            return;
        end
        % Find Subject directories in selected path
        dir_path = dir(selpath);
        dir_path(contains({dir_path.name},{'output'})) = [];
        dir_path(contains({dir_path.name},{'.'})) = [];
        dir_path = {dir_path.name};
    else
        [fn,shortpath,flag,LOC_turbo_path] = checkTurboPath(varargin{1});
        if all(flag)
            warning('Selected path is not on Turbo storage.');
            return;
        end
    end
    nf = numel(dir_path);
    
    % Determine SBATCH inputs
    tstr = char(datetime('now','Format','yyyyMMddHHmmss'));
    jobname = sprintf('GL_vesselSeg_%s',tstr);
    fname = [jobname,'.sh'];
    part_str = 'standard';
    pmem = 16; % max GB per process
    ptime = 60; % max minutes per process
    mxmem = 180; % 180GB max memory for standard node
    cores = min(min(nf,floor(mxmem/pmem))+1,36);
    walltime = ptime*ceil(nf/(cores-1));
    if N==1
        Nnodes = ceil(walltime/(14*24*60)); % 14 days max walltime
    else
        Nnodes = N;
    end
    % Recalculate taking into account #nodes
    nf = ceil(nf/Nnodes);
    cores = min(min(nf,floor(mxmem/pmem))+1,36);
    mem = pmem * (cores-1);
    walltime = ptime*ceil(nf/(cores-1));
    if walltime > (14*24*60) % 14 day max walltime
        warning('Processes may not complete.\nEstimated time = %s\n',...
            datestr(duration(0,walltime,0),'dd-HH:MM:SS'));
    end
    
    % Determine batch job inputs and save input file for GL execution:
    tempdir = fullfile(LOC_turbo_path,'GreatLakes','temp');
    fname_inputs = fullfile(tempdir,sprintf('batch_inputs_%s.mat',tstr));
    GL_turbo_dir = {'/nfs/turbo/umms-cgalban'};
    dir_path = cellfun(@(x,y)strjoin([GL_turbo_dir,shortpath(2:end),fn,{y}],'/'),dir_path);
    save(fname_inputs,'dir_path','Nnodes');
    
    % Write SBATCH file
    opts = {'partition',part_str, 'cores',cores, 'mem',mem, 'walltime',walltime,...
        'username',username, 'inputs_fname',fname_inputs};
    if Nnodes > 1
        opts = [opts,{'array',Nnodes}];
    end
    write_sbatch(fullfile(tempdir,fname),jobname,GLstr,opts{:});
    
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
    
    if ischar(varargin{2}) && exist(varargin{2},'file') && strcmp(varargin{2}(end-3:end),'.mat')
        inputs_fname = varargin{2};
    else
        error('Invalid input: must be valid .mat file containing function inputs.');
    end
    
    fprintf('Available workers on this node: %s\n',getenv('SLURM_CPUS_PER_TASK'));
    
    p = load(inputs_fname);
    flds = {'dir_path'};
    ind = ~isfield(p,flds);
    if any(ind)
        error(['Invalid input. Missing:',repmat(' %s',1,nnz(ind))],flds{ind});
    elseif ~(iscellstr(p.fn_ins) && iscellstr(p.fn_seg) && ischar(p.sv_path) ...
            && (numel(p.fn_ins)==numel(p.fn_seg)))
        error('Invalid input. fn_ins and fn_seg must be cellstr with equal number of elements.');
    end
    if isfield(p,'Nnodes')
        njobs = p.Nnodes;
    else
        njobs = 1;
    end
    
    % Check for job array:
    slurm = getSLURM;
    jobstr = sprintf('%s_array%u',slurm.jobname,slurm.arraytaskid);
    if ~isempty(slurm.arraytaskid)
        fprintf('Running job number %u of %u\n',slurm.arraytaskid,njobs);
    end
    nf = numel(p.dir_path);
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
    
    pp = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(pp)
        parpool(4);
    else
        poolsize = pp.NumWorkers
    end
    pp = gcp();
    
    dir_path = p.dir_path;
    
    T = [];
    parfor procdir = string(dir_path)
        tT = tPRM_aim2_sub(procdir);
        T = [T;tT];
    end
    save(fullfile(fileparts(dir_path{1}),sprintf('tPRM_aim2_Results_%s',jobstr)),'T');
end