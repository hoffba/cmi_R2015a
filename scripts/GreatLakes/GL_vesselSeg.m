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
        pmem = 12; % max GB per process
        ptime = 300; % max minutes per process
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
        jobname = getenv('SLURM_JOB_NAME');
        jobnum = str2double(getenv('SLURM_ARRAY_TASK_ID'));
        if ~isnan(jobnum)
            fprintf('Running job number %u of %u\n',jobnum,njobs);
        end
        nf = numel(p.fn_ins);
        if ~isnan(jobnum)
            ind = round((nf/njobs)*[jobnum-1 jobnum]);
            
            fprintf('ind: %u-%u\n',ind+[1,0]);
            
            ind = (ind(1)+1):ind(2);
        else
            ind = 1:nf;
        end
        nf = numel(ind);
        
        % Set up cluster properties
        c = parcluster;
        jobdir = fullfile(c.JobStorageLocation,sprintf('%s_array%u',jobname,jobnum));
        mkdir(jobdir);
        c.JobStorageLocation = jobdir;
        
        % Start batch jobs
        fn_base = cell(1,nf);
        for i = 1:nf
            % Find subsection index:
            ii = ind(i);
            
            % Extract base name
            tok = regexp(p.fn_ins{ii},'/([^./]+)\.','tokens');
            fn_base{i} = tok{1}{1};
            
            fprintf('Starting batch job %u:\n',i,p.fn_ins{ii});
            fprintf('   Ins: %s\n   Seg: %s\n   Sv: %s\n',p.fn_ins{ii},p.fn_seg{ii},p.sv_path);
            job(i) = batch(c,@vesselSeg_BH,1,[p.fn_ins(ii),p.fn_seg(ii),p.sv_path]);
        end
                
        % Wait for all jobs to complete
        errflag = true(nf,1);
        dt = zeros(nf,1);
        for i = 1:nf
            wait(job(i));
            dt(i) = minutes(job(i).FinishDateTime - job(i).StartDateTime);
            fprintf('Job %u finished after %.1f minutes.\n',i,dt(i));
            if strcmp(job(i).State,'failed')
                errflag(i) = false;
                fprintf(job(i).Tasks(1).ErrorMessage);
            end
                        
            % Delete job files:
            job(i).delete;
            
        end
        rmdir(jobdir,'s');
        
        fprintf('Processing complete.\nAverage processing time = %.1f (%.1f) minutes.\n',...
            mean(dt(errflag)),std(dt(errflag)));
        
    else
        error('Invalid input');
    end

end

function val = job_vesselSeg(subj_dir)
% Perform lung vessel segmentation 
% Inputs: subj_dir = directory containing all subject data

    [~,ID] = fileparts(subj_dir);
    
    % Find and load INSP CT file:
    fn_ins = dir(fullfile(subj_dir,'*_INSP_*.ct.nii.gz'));
    if isempty(fn_ins), warning('No INSP file found for subject: %s.',ID); return, end
    fn_ins = {fn_ins.name};
    ind = startsWith(fn_ins,'re_');
    re_flag = false;
    if any(ind)
        fn_ins = fn_ins{ind};
    else
        fn_ins = fn_ins{1};
        re_flag = true;
    end
    info_ct = niftiinfo(fn_ins);
    ct = double(niftiread(info_ct));
    if re_flag
        [ct,info_ct] = resample_subj(ct,info_ct,fn_ins);
    end
    
    % Find INSP segmentation file:
    fn_seg = dir(fullfile(subj_dir,'*_INSP_*.lobe_segmentation.nii.gz'));
    if isempty(fn_seg), warning('No INSP segmentation file found for subject: %s.',ID); return, end
    fn_seg = {fn_seg.name};
    ind = startsWith(fn_seg,'re_');
    re_flag = false;
    if any(ind)
        fn_seg = fn_seg{ind};
    else
        fn_seg = fn_seg{1};
        re_flag = true;
    end
    info_seg = niftiinfo(fn_seg);
    seg = niftiread(info_seg);
    if re_flag
        [seg,info_seg] = resample_subj(seg,info_seg,fn_seg);
    end
    
    % Save eroded lobe map:
    lobe_tag = [11,12,13,21,22];
    se = strel('sphere',5);
    eroded_lobes = zeros(size(seg));
    for i = 1:numel(lobe_tag)
        eroded_lobes(imerode(seg == lobe_tag(i),se)) = lobe_tag(i);
    end
    eroded_lobes = uint8(eroded_lobes);
    save(fullfile(subj_dir,[ID,'_erodedLobes.mat']),'eroded_lobes');
    eroded_lobes = eroded_lobes > 0;
    
    % Full lung volume:
    segBW = ismember(seg,[11,12,13,21,22]);
    
    % Determine temp folder for fast saving
%     fn_base = regexp(fn_ins,'/([^./]+)\.','tokens');
%     fn_base = fn_base{1}{1};
%     if isfolder('/tmpssd')
%         tmpdir = '/tmpssd';
%     elseif isfolder('/tmp')
%         tmpdir = '/tmp';
%     else
%         tmpdir = sv_path;
%     end

%     % Set up labels for results
%     info_sv = info_ins;
%     info_sv.Datatype = 'int16';
%     info_sv.BitsPerPixel = 16;

    % Generate ero

    % Generate enhanced vessel maps
    fprintf('Calculating enhanced vessel maps ...\n');
    [vessels, ~] = MTHT3D( Normalize(ct.*segBW) ,...
                     0.5:1:5.5 ,...
                     12 ,... % number of orientations
                     70 ,... % beta
                     0.5 ,... % alpha
                     15 ,... % c : parameters for Vesselness
                     -1/3 );  % alfa : parameters for Neuriteness
    vessels = double(vessels) .* segBW;
    save(fullfile(subj_dir,[ID,'_enhancedVessel.mat']),'vessels');

    % Binarize vessels:
    info = info_seg;
    info.Datatype = 'int8';
    info.BitsPerPixel = 8;
    bin_vessels = activecontour(vessels,imbinarize(vessels,'adaptive').*eroded_lobes,5,'Chan-Vese');
    niftiwrite(uint8(bin_vessels),fullfile(subj_dir,[ID,'_binVessel.nii']),'Compressed',true);
    
    % Generate CSA maps
    csa_map = createCSA_maps(bin_vessels);
    save(fullfile(subj_dir,[ID,'_CSA_skel.mat']),'csa_map');

    %Frangi Filter
    opt = struct('FrangiScaleRange',[0.5 5.5],...
                 'FrangiScaleRatio',1,...
                 'BlackWhite',false);
    frangi_enhanced_vessels = FrangiFilter3D(ct.*segBW,opt); 
    save(fullfile(subj_dir,[ID,'_enhancedVessel_frangi.mat']),'frangi_enhanced_vessels');

    %Curvilinear Filter
    curv_enhanced_vessels = vesselness3D(ct.*segBW, 0.5:1:5.5, [1,1,1], 1, true).*segBW;
    save(fullfile(subj_dir,[ID,'_enhancedVessel_curv.mat']),'curv_enhanced_vessels');
    
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

function [I,info] = resample_subj(I,info,fname)
    
    if endsWith(fname,'.lobe_segmentation.nii.gz')
        interpmethod = 'nearest';
    else
        interpmethod = 'linear';
    end

    %Get new size
    pixDim = info.PixelDimensions;
    pixSize = 0.625;
    
    [sx,sy,sz] = size(I);
    xq = length(0:pixSize/pixDim(1):sx);
    yq = length(0:pixSize/pixDim(2):sy);
    zq = sz; %was length(0:sz) so original batch has one extra z layer
    
    I = imresize3(img_ct, [xq yq zq], 'method', interpmethod);
    
    info.PixelDimensions = [0.625, 0.625, 0.625];
    info.ImageSize = size(I);
    
    % Save resampled data to file:
    [subj_dir,fname,ext] = fileparts(fname);
    fname = fullfile(subj_dir,['re_',fname,ext]);
    niftiwrite(I,fname,info);
    
end
