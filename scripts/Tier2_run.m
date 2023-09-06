function Tier2_run(varargin)
% Prepares processes for running on the Tier2 server: galban-ap-ps1a
% Saves temporary files to: Turbo:\Tier2\temp
% Must execute shell script manually using PuTTy
% Syntax:
%   GL_run( username , ...   % your uniquename (e.g. 'cgalban')
%           function_string , ... % name of recursive function (e.g. 'job_tPRM')
%           inputs , ...     % <1xM> cell array - inputs to recursive function
%           path_flag , ...  % <1xM> logical - (optional) flag to correct Turbo path (local to GL)
%           static_flag , ...% <1xM> logical - indicateds static input accross increments
%           Name/Value )
% Optional Name/Value inputs:
%   'Nodes' - number of nodes to run on
%   'ProcessTime' - estimated processing time (minutes) per iteration (overestimate)
%   'ProcessMemory' - estimated memory requirement (GB) per iteration
%   'Partition' - type of partition to use {'auto','standard','largemem'}
% Example:
%   GL_run('bahoff','job_tPRM',{fn_prm,fn_seg,sv_path},[true,true,true],[false,false,true],...
%          'Nodes',5,'ProcessTime',120,'ProcessMemory',12,'Partition','auto')

if nargin~=2
    
% ~~~~~~~~~~~~~~~~~~~~~
% ~~ LOCAL EXECUTION ~~
% ~~~~~~~~~~~~~~~~~~~~~
    
    opts = parseInputs(varargin{:});

    % Determine SBATCH inputs
    jobname = sprintf('%s_%s_%s',opts.username,opts.function_string,opts.TimeStamp);
    fname = [jobname,'.sh'];

    % Save batch job inputs to Turbo temp directory
    fname_inputs = fullfile(opts.TempDir,sprintf('%s_INPUTS.mat',jobname));
    save(fname_inputs,'-struct','opts');

    % Write shell file
    sh_opts = {'username',opts.username, 'inputs_fname',fname_inputs};
    write_shellscript(fullfile(opts.TempDir,fname),jobname,'Tier2_run',sh_opts{:});

    % Copy command to run on Tier2 server terminal
    sb_str = sprintf('cd /nfs/turbo/umms-cgalban/Tier2/temp && sbatch %s',fname);
    clipboard('copy',sb_str);
    fprintf(['\n\n1) Use PuTTy to log onto Tier2:\n',...
             '   Host Name: galban-ap-ps1a\n',...
             '   Login: Level 2\n',...
             '2) Start shell script using terminal commands:\n',...
             '   (ALREADY IN CLIPBOARD, just paste into PuTTy)\n',...
             '   ',sb_str,'\n'],fname);

else
    
% ~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~ TIER2 EXECUTION ~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~
    
        % In Tier2, load file containing relevant inputs to the processing function
        %   and start processing jobs  
        
        if ischar(varargin{2}) && exist(varargin{2},'file') && strcmp(varargin{2}(end-3:end),'.mat')
            inputs_fname = varargin{2};
        else
            error('Invalid input: must be valid .mat file containing function inputs.');
        end
        
        fprintf('Available workers on this node: %s\n',getenv('SLURM_CPUS_PER_TASK'));
        
        p = load(inputs_fname);
        
        % Check for job array:
        jobname = getenv('SLURM_JOB_NAME');
        jobnum = str2double(getenv('SLURM_ARRAY_TASK_ID'));
        if isnan(jobnum)
            ind = 1:p.Niter;
        else
            fprintf('Running job number %u of %u\n',jobnum,p.Nodes);
            
            ind = round(p.Niter/p.Nodes * [jobnum-1 jobnum]) + [1,0];
            fprintf('ind: %u-%u\n',ind);
            
            ind = ind(1):ind(2);
        end
        p.Niter = numel(ind);
        
        % Set up cluster properties
        c = parcluster;
        jobdir = fullfile(c.JobStorageLocation,sprintf('%s_array%u',jobname,jobnum));
        mkdir(jobdir);
        c.JobStorageLocation = jobdir;
        
        % Start batch jobs
        Nin = numel(p.inputs);
        for i = 1:p.Niter
            % Find subsection index:
            ii = ind(i);
            
            % Extract inputs:
            temp_inputs = cell(1,Nin);
            for j = 1:Nin
                if p.static_flag(j)
                    temp_inputs{j} = p.inputs{j};
                elseif iscell(p.inputs{j})
                    temp_inputs{j} = p.inputs{j}{ii};
                else
                    temp_inputs{j} = p.inputs{j}(ii);
                end
            end
            
            fprintf('Starting batch job %u\n',ii);
            job(i) = batch(c,str2func(p.function_string),1,temp_inputs);
        end
        
        % Initialize table for statistics:
        T = table;
        
        % Wait for all jobs to complete
        errflag = true(p.Niter,1);
        dt = zeros(p.Niter,1);
        for i = 1:p.Niter
            wait(job(i));
            if ~isempty(job(i).Tasks(1).Error) || strcmp(job(i).State,'failed')
                errflag(i) = false;
                try
                    getReport(job(i).Tasks(1).Error)
                catch
                    job(i).Tasks(1).Error
                end
            else
                % Compile statistics:
                dt(i) = minutes(job(i).FinishDateTime - job(i).StartDateTime);
                val = job(i).fetchOutputs;
                val = val{1};
                if istable(val) && ~isempty(val) && (isempty(T) || (size(val,2) == size(T,2)))
                    val.Properties.RowNames = {};
                    T = [T;val];
                elseif ischar(val)
                    fprintf(val);
                end
            end
            fprintf('Job %u finished after %.1f minutes.\n',i,dt(i));
            fprintf('  State: %s\n',job(i).State);
            job(i).diary
            
            % Delete job files:
            job(i).delete;
            
        end
        rmdir(jobdir,'s');
        
        % Write table with mean values to save path
        if ~isempty(p.save_path)
            svdir = checkTurboPath(p.save_path);
            jobnum_str = '';
            if ~isnan(jobnum)
                jobnum_str = [num2str(jobnum),'_'];
            end
            svname = fullfile(svdir{1},sprintf('%s_%sResults',jobname,jobnum_str));
            fprintf('Attempting to save results:\n');
            if istable(T)
                svname = [svname,'.csv'];
                writetable(T,svname,'WriteRowNames',true);
            else
                svname = [svname,'.mat'];
                save(svname,T);
            end
            fprintf('   ... %s\n',svname);
        end
        
        fprintf('Processing complete.\nAverage processing time = %.1f (%.1f) minutes.\n',...
            mean(dt(errflag)),std(dt(errflag)));
end


function opts = parseInputs(varargin)
p = inputParser;
addRequired(p,'username',@ischar);
addRequired(p,'function_string',@(x)ischar(x)&&(exist(x,'file')==2));
addRequired(p,'inputs',@iscell);
addRequired(p,'path_flag',@islogical);
addRequired(p,'static_flag',@islogical);
addParameter(p,'Nodes',1,@isscalar);
addParameter(p,'ProcessTime',60,@isscalar);
addParameter(p,'ProcessMemory',6,@isscalar);
addParameter(p,'Partition','auto',@(x)ismember(x,{'auto','standard','largemem'}));
addParameter(p,'mods',{},@iscellstr);
addParameter(p,'TimeStamp','',@ischar);
addParameter(p,'save_path','',@ischar);
parse(p,varargin{:});
opts = p.Results;

% Validate recursive inputs
opts.Niter = numel(opts.inputs{find(~opts.static_flag,1)});
for i = 1:numel(opts.inputs)
    % Check agreement in number of iterations
    if ~opts.static_flag(i) && numel(opts.inputs{i})~=opts.Niter
        error('Number of iterations must agree between inputs.');
    end
    
    % Validate directory inputs:
    if opts.path_flag(i)
        tclass = class(opts.inputs{i});
        [opts.inputs{i},~,LOC_Turbo_Path] = checkTurboPath(opts.inputs{i});
        if ismember(tclass,{'char','string'})
            opts.inputs{i} = opts.inputs{i}{1};
        end
    end
end

% Determine what type of partition to use
partition_stats = {'standard',180;...
                   'largemem',1500};
if strcmp(opts.Partition,'auto')
    if opts.ProcessMemory > 20
        opts.Partition = 'largemem';
    else
        opts.Partition = 'standard';
    end
end
opts.MaxMem = partition_stats{strcmp(opts.Partition,partition_stats(:,1)),2};

opts.TempDir = fullfile(LOC_Turbo_Path,'GreatLakes','temp');


function write_shellscript()
