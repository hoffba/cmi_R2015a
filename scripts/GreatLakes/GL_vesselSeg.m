function [fn_ins,fn_seg,sv_path] = GL_vesselSeg(varargin)
    % Performs tPRM on existing PRM maps on GreatLakes
    % Syntax:
    %   GL_tPRM(username)
    %   GL_tPRM(username,N)
    %   GL_tPRM(username,fn_ins,fn_seg,sv_path)
    %   GL_tPRM(username,fn_ins,fn_seg,sv_path,N)
    % Required inputs:
    %       fn_ins = cell array of Inspiratory CT files (*.nii.gz)
    %       fn_seg = cell array of Inspiratory segmentation files (*.nii.gz)
    %             OR string location of folder containing matching seg files
    %       sv_path = string location of directory for saving results
    %       N = requested node array

    % Parse inputs
    if ischar(varargin{1}) % Should be uniquename
        
% ~~~~~~~~~~~~~~~~~~~~~
% ~~ LOCAL EXECUTION ~~
% ~~~~~~~~~~~~~~~~~~~~~
        
        % Parse Inputs:
        N = 1;
        username = varargin{1};
        gflag = false; % GUI flag
        if nargin>1
            if isnumeric(varargin{2}) && isscalar(varargin{2})
                N = varargin{2};
                gflag = true;
            elseif (nargin>3) && (ischar(varargin{2}) || iscellstr(varargin{2}))
                fn_ins = varargin{2};
                fn_seg = varargin{3};
                sv_path = varargin{4};
                if (nargin==5) && isnumeric(varargin{5}) && isscalar(varargin{5})
                    N = varargin{5};
                end
            else
                warning('Invalid input #2.');
                return;
            end
        else
            gflag = true;
        end
            
        % Select data for processing:
        if gflag
            % User GUI for file selection:
            [fn_ins,fpath_ins] = uigetfile('*_INSP_*.nii.gz','Select INS data for processing:','MultiSelect','on');
            if ~fn_ins,return;end
            fn_ins = fullfile(fpath_ins,fn_ins);
            
            fn_seg = uigetdir(pwd,'Select folder containing corresponding lobe_segmentation.nii.gz files:');
            if ~fn_seg,return;end
            
            sv_path = uigetdir(pwd,'Select folder for saving vessel results:');
            if ~sv_path,return;end
        else
            % Validate fn inputs
            if ~(iscellstr(fn_ins) || ischar(fn_ins))
                error('Input 2 (fn_ins) must be a string or cellstr containing INSP file locations.');
            end
            if ~(iscellstr(fn_seg) || ischar(fn_seg))
                error('Input 3 (fn_seg) must be a directory containing segmentation files or cellstr of individual file locations.');
            end
            if ~ischar(sv_path)
                error('Input 4 (sv_path) must be a directory location for saving results.');
            end
        end
        if ischar(fn_ins)
            fn_ins = {fn_ins};
        end
        folder_flag = false;
        if ischar(fn_seg)
            folder_flag = isfolder(fn_seg);
            fn_seg = {fn_seg};
        end
        if ~folder_flag && (numel(fn_ins)~=numel(fn_seg))
            error('Each INSP file must have one corresponding segmentation file.');
        end
        
        % Validate save directory:
        netdrives = findNetDrives('umms-cgalban');
        LOC_turbo_path = netdrives.Drive;
        [~,shortpath_sv,flag] = checkTurboPath(sv_path,LOC_turbo_path);
        if flag
            error('Save directory must be on Turbo.');
        end
        
        % Validate INSP files:
        [fn_ins,shortpath_ins,flag_ins] = checkTurboPath(fn_ins,LOC_turbo_path);
        if all(flag_ins)
            error('No valid INSP files were found.');
        elseif ~folder_flag
            fn_seg(flag_ins) = [];
        end
        
        % Validate segmentation files:
        [fn_seg,shortpath_seg,flag_seg] = checkTurboPath(fn_seg,LOC_turbo_path);
        if all(flag_seg)
            error('No valid segmentation files were found.');
        elseif ~folder_flag
            if any(flag_seg)
                % Show files that weren't found
                fprintf('%u segmentation files not found:\n',nnz(flag_seg));
                fprintf('   %s\n',fn_seg{flag_seg});
                fn_ins(flag_seg) = [];
            end
        end
        
        % Find corresponding segmentation files
        if folder_flag
            nf = numel(fn_ins);
            fn_seg = cell(nf,1);
            flag_seg = false(nf,1);
            for i = 1:nf
                tname = dir(fullfile(shortpath_seg{1}{:},...
                    [extractBefore(fn_ins{i},'.'),'*.lobe_segmentation.nii.gz']));
                if isempty(tname)
                    flag_seg(i) = true;
                else
                    fn_seg{i} = tname.name;
                end
            end
            if all(flag_seg)
                error('No valid segmentation files were found.');
            elseif any(flag_seg)
                % Show files that weren't found
                fprintf('%u segmentation files not found:\n',nnz(flag_seg));
                fprintf('   %s\n',fn_seg{flag_seg});
                fn_seg(flag_seg) = [];
            end
        end
        
        % Determine SBATCH inputs
        tstr = char(datetime('now','Format','yyyyMMddHHmmss'));
        jobname = sprintf('GL_vesselSeg_%s',tstr);
        fname = [jobname,'.sh'];
        part_str = 'standard';
        pmem = 32; % max GB per process
        ptime = 720; % max minutes per process
        mxmem = 180; % 180GB max memory for standard node
        nf = numel(fn_ins);
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
        fn_ins = cellfun(@(x,y)strjoin([GL_turbo_dir,x(2:end),{y}],'/'),...
            shortpath_ins,fn_ins,'UniformOutput',false);
        fn_seg = cellfun(@(x,y)strjoin([GL_turbo_dir,x(2:end),{y}],'/'),...
            shortpath_seg,fn_seg,'UniformOutput',false);
        sv_path = strjoin([GL_turbo_dir,shortpath_sv{1}(2:end)],'/');
        save(fname_inputs,'fn_ins','fn_seg','sv_path','Nnodes');
        
        % Write SBATCH file
        opts = {'partition',part_str, 'cores',cores, 'mem',mem, 'walltime',walltime,...
            'username',username, 'inputs_fname',fname_inputs};
        if Nnodes > 1
            opts = [opts,{'array',Nnodes}];
        end
        write_sbatch(fullfile(tempdir,fname),jobname,'GL_vesselSeg',opts{:});
        
        % Run SBATCH
        sb_str = sprintf('cd /nfs/turbo/umms-cgalban/GreatLakes/temp && sbatch %s',fname);
        clipboard('copy',sb_str);
        fprintf(['\n\n1) Use PuTTy to log onto Great Lakes:\n',...
                 '   Host Name: greatlakes.arc-ts.umich.edu\n',...
                 '   Login: Level 1\n',...
                 '2) Start SBATCH using terminal commands:\n',...
                 '   (ALREADY IN CLIPBOARD, just paste into PuTTy using middle mouse button)\n',...
                 '   ',sb_str,'\n'],fname);
        
    elseif (nargin==2) % ( 1 , inputs_filename )
        
% ~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~ GREATLAKES EXECUTION ~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        % In GL, load file containing relevant inputs to the processing function
        %   and start processing jobs  
        
        if ischar(varargin{2}) && exist(varargin{2},'file') && strcmp(varargin{2}(end-3:end),'.mat')
            inputs_fname = varargin{2};
        else
            error('Invalid input: must be valid .mat file containing function inputs.');
        end
        
        fprintf('Available workers on this node: %s\n',getenv('SLURM_CPUS_PER_TASK'));
        
        p = load(inputs_fname);
        flds = {'fn_ins','fn_seg','sv_path'};
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
        if ~isempty(slurm.arraytaskid)
            fprintf('Running job number %u of %u\n',slurm.arraytaskid,njobs);
        end
        nf = numel(p.fn_ins);
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
        jobdir = fullfile(c.JobStorageLocation,sprintf('%s_array%u',slurm.jobname,slurm.arraytaskid));
        if ~isfolder(jobdir)
            mkdir(jobdir);
        end
        c.JobStorageLocation = jobdir;
        
        % Start batch jobs
        ID = cell(nf,1);
        for i = 1:nf
            % Find subsection index:
            ii = ind(i);
            
            [~,ID{i}] = fileparts(p.fn_ins{ii});
            if startsWith(ID{i},'re_')
                ID{i}(1:3) = [];
            end
            if contains(ID{i},'_')
                ID{i} = extractBefore(ID{i},'_');
            elseif contains(ID,'.')
                ID{i} = extractBefore(ID{i},'.');
            end

            % Save directly to turbo
            svdir = fullfile(p.sv_path,ID{i});
            
            fprintf('Starting batch job %u:\n',i);
            fprintf('   Ins: %s\n   Seg: %s\n   Sv: %s\n',p.fn_ins{ii},p.fn_seg{ii},p.sv_path);
            job(i) = batch(c,@vesselSeg_BH,1,[p.fn_ins(ii),p.fn_seg(ii),svdir]);
        end
                
        % Wait for all jobs to complete
        errflag = true(nf,1);
        dt = zeros(nf,1);
        T = [];
        for i = 1:nf
            wait(job(i));
            
            dt(i) = minutes(datetime('now','TimeZone','local') - job(i).StartDateTime);
            
            if isempty(job(i).Tasks(1).Error)
                fprintf('Job %u finished after %.1f minutes.\n',i,dt(i));
                tT = fetchOutputs(job(i));
                T = [T;tT{1}];
            else
                fprintf('Job %u failed after %.1f minutes.\n',i,dt(i));
                errflag(i) = false;
                fprintf('%s\n',job(i).Tasks(1).ErrorMessage);
                job(i).Tasks(1)
                fprintf('File: %s\n',job(i).Tasks(1).stack(1).file);
                fprintf('Name: %s\n',job(i).Tasks(1).stack(1).name);
                fprintf('Line: %s\n',job(i).Tasks(1).stack(1).line);
            end
                    
            job(i).diary
            
            % Delete job files:
            job(i).delete;
            
        end
        rmdir(jobdir,'s');
        
        fprintf('Processing complete.\nAverage processing time = %.1f (%.1f) minutes.\n',...
            mean(dt(errflag)),std(dt(errflag)));
        
        tfname = fullfile(p.sv_path,sprintf('vesselSeg_Results_%s.csv',slurm.jobname));
        fprintf('Saving tabulated results: %s',tfname);
        writetable(T,tfname);
        
    else
        error('Invalid input');
    end

end

function [fn,shortpath,flag] = checkTurboPath(fnames,dcode)
    if ischar(fnames)
        fnames = {fnames};
    end
    nf = numel(fnames);
    flag = false(nf,1); % flag fnames to remove
    fn = cell(nf,1);
    shortpath = cell(nf,1);
    for i = 1:nf
        str = strsplit(fnames{i},'\');
        if ~strcmp(str{1},dcode)
            flag(i) = true;
            fprintf('File not on Turbo: %s\n',fnames{i});
        else
            if isfolder(fnames{i})
                if nf == 1
                    shortpath{i} = str;
                else
                    flag(i) = true;
                    fprintf('Value must be a file path, not directory: %s\n',fnames{i});
                end
            else
                shortpath{i} = str(1:end-1);
                fn(i) = str(end);
            end
        end
    end
    fn(flag) = [];
    shortpath(flag) = [];
end

