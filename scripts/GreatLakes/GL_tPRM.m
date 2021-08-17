function [fn_prm,fn_seg,sv_path] = GL_tPRM(varargin)
    % Performs tPRM on existing PRM maps on GreatLakes
    % Syntax:
    %   GL_tPRM(username)
    %   GL_tPRM(username,fn_prm,fn_seg,sv_path)
    % Required inputs:
    %       fn_prm = cell array of PRM files (*.nii.gz)
    %       fn_seg = cell array of segmentation files (*.nii.gz)
    %             OR string location of folder containing matching seg files
    %       sv_path = string location of directory for saving results

    % Parse inputs
    switch nargin
        case 1
            if isfile(varargin{1})
                % GL EXECUTION
                flag = false;
                inputs_fname = varargin{1};
            else
                flag = true;
                username = varargin{1};
            end
        case 4
            flag = true;
            username = varargin{1};
            fn_prm = varargin{2};
            fn_seg = varargin{3};
            sv_path = varargin{4};
        otherwise
            error('Invalid number of inputs.');
    end
    
    if flag
        
        % ~~ LOCAL EXECUTION ~~
        
        % Select data for processing:
        if nargin==4
            if ~(iscellstr(fn_prm) || ischar(fn_prm))
                error('Input 1 (fn_prm) must be a string or cellstr containing PRM file locations.');
            end
            if ~(iscellstr(fn_seg) || ischar(fn_seg))
                error('Input 2 (fn_seg) must be a directory containing segmentation files or cellstr of individual file locations.');
            end
            if ~ischar(sv_path)
                error('Input 3 (sv_path) must be a directory location for saving results.');
            end
        else
            % User GUI for file selection:
            [fn_prm,fpath_prm] = uigetfile('*.prm.nii.gz','Select PRM data for processing:','MultiSelect','on');
            if isempty(fn_prm)
                return;
            end
            fn_prm = fullfile(fpath_prm,fn_prm);
            fn_seg = uigetdir(pwd,'Select folder containing corresponding lobe_segmentation.nii.gz files:');
            sv_path = uigetdir(pwd,'Select folder for saving tPRM results:');
        end
        if ischar(fn_prm)
            fn_prm = {fn_prm};
        end
        folder_flag = false;
        if ischar(fn_seg)
            folder_flag = isfolder(fn_seg);
            fn_seg = {fn_seg};
        end
        if ~folder_flag && (numel(fn_prm)~=numel(fn_seg))
            error('Each PRM file must have one corresponding segmentation file.');
        end
        
        % Validate save directory:
        netdrives = findNetDrives('umms-cgalban');
        LOC_turbo_path = netdrives.Drive;
        [~,shortpath_sv,flag] = checkTurboPath(sv_path,LOC_turbo_path);
        if flag
            error('Save directory must be on Turbo.');
        end
        
        % Validate PRM files:
        [fn_prm,shortpath_prm,flag_prm] = checkTurboPath(fn_prm,LOC_turbo_path);
        if all(flag_prm)
            error('No valid PRM files were found.');
        elseif ~folder_flag
            fn_seg(flag_prm) = [];
        end
        
        % Validate segmentation files:
        [fn_seg,shortpath_seg,flag_seg] = checkTurboPath(fn_seg,LOC_turbo_path);
        if all(flag_seg)
            error('No valid segmentation files were found.');
        elseif ~folder_flag
            fn_prm(flag_seg) = [];
        end
        
        % Find corresponding segmentation files
        if folder_flag
            fn_seg = cellfun(@(x)regexprep(x,'.prm.','.njh.lobe_segmentation.'),fn_prm,'UniformOutput',false);
            fnames = fullfile(shortpath_seg{1}{:},fn_seg);
            ind = ~cellfun(@(x)exist(x,'file'),fnames);
            if all(ind)
                error('No valid segmentation files were found.');
            elseif any(ind)
                % Show files that weren't found
                fprintf('%u segmentation files not found:\n',nnz(ind));
                fprintf('   %s\n',fn_seg{ind});
                fn_seg(ind) = [];
            end
        end

        % Determine batch job inputs and save input file for GL execution:
        tstr = char(datetime('now','Format','yyyyMMddHHmmss'));
        tempdir = fullfile(LOC_turbo_path,'GreatLakes','temp');
        fname_inputs = sprintf('batch_inputs_%s.mat',tstr);
        GL_turbo_dir = {'/nfs/turbo/umms-cgalban'};
        fn_prm = cellfun(@(x,y)strjoin([GL_turbo_dir,x(2:end),{y}],'/'),...
            shortpath_prm,fn_prm,'UniformOutput',false);
        fn_seg = cellfun(@(x,y)strjoin([GL_turbo_dir,x(2:end),{y}],'/'),...
            shortpath_seg,fn_seg,'UniformOutput',false);
        sv_path = strjoin([GL_turbo_dir,shortpath_sv{1}(2:end)],'/');
        save(fullfile(tempdir,fname_inputs),'fn_prm','fn_seg','sv_path');
        
        % Write SBATCH file
        nf = numel(fn_prm);
        fname = sprintf('GL_tPRM_%s.sh',tstr);
        cores = min(nf,6);
        mem = cores * 24;
        write_sbatch(fullfile(tempdir,fname),fname,username,cores,mem,'GL_tPRM',...
            strjoin([GL_turbo_dir,{'GreatLakes','temp',fname_inputs}],'/'));
        
        % Run SBATCH
        clipboard('copy',sprintf('cd /nfs/turbo/umms-cgalban/GreatLakes/temp && sbatch %s',fname));
        fprintf(['\n\n1) Use PuTTy to log onto Great Lakes:\n',...
                 '   Host Name: greatlakes.arc-ts.umich.edu\n',...
                 '   Login: Level 1\n',...
                 '2) Start SBATCH using terminal commands:\n',...
                 '   (ALREADY IN CLIPBOARD, just paste into PuTTy using middle mouse button)\n',...
                 '   cd /nfs/turbo/umms-cgalban/GreatLakes/temp\n',...
                 '   sbatch %s\n'],fname);
        
    elseif ischar(inputs_fname) && exist(inputs_fname,'file') && strcmp(inputs_fname(end-3:end),'.mat')
        
        % ~~ GREATLAKES EXECUTION ~~
        
        % In GL, load file containing relevant inputs to the processing function
        %   and start processing jobs  
        p = load(inputs_fname);
        flds = {'fn_prm','fn_seg','sv_path'};
        ind = ~isfield(p,flds);
        if any(ind)
            error(['Invalid input. Missing:',repmat(' %s',1,nnz(ind))],flds{ind});
        elseif ~(iscellstr(p.fn_prm) && iscellstr(p.fn_seg) && ischar(p.sv_path) ...
                && (numel(p.fn_prm)==numel(p.fn_seg)))
            error('Invalid input. fn_prm and fn_seg must be cellstr with equal number of elements.');
        end
        
        % Start batch jobs
        nf = length(p.fn_prm);
        for i = 1:nf
            fprintf('Starting batch job %u:\n   %s\n',i,p.fn_prm{i});
            job(i) = batch(@job_tPRM,0,[p.fn_prm(i),p.fn_seg(i),p.sv_path]);
        end
        
        % Wait for all jobs to complete
        for i = 1:nf
            wait(job(i),'finished');
            fprintf('Job %u finished.\n',i);
            diary(job(i))
%             getReport(job(i).Tasks(1).Error)
        end
        
    else
        error('Invalid input');
    end

end

function job_tPRM(fn_prm,fn_seg,sv_path)
% Perform tPRM processing 
% Inputs: fn_prm = full filename of PRM classification map (as seen by GL)
%         fn_seg = full filename of segmentation map (as seen by GL)
%         sv_path = path to save resulting tPRM maps into
    
    fn_base = regexp(fn_prm,'/([^./]+)\.','tokens');
    fn_base = fn_base{1}{1};

    % Load PRM map
    fprintf('Loading PRM map: %s\n',fn_prm);
    info_prm = niftiinfo(fn_prm);
    prm = niftiread(info_prm);
    voxsz = info_prm.PixelDimensions;

    % Load lobe segmentation
    fprintf('Loading segmentation: %s\n',fn_seg);
    info_seg = niftiinfo(fn_seg);
    seg = niftiread(info_seg);
    % Only need total lung segmentation for this:
    seg = ismember(seg,[11,12,13,21,22]);

    % Set up labels for results
    info_sv = info_prm;
    info_sv.Datatype = 'int16';
    info_sv.BitsPerPixel = 16;
    prmclass = struct('val',{3,(4:5),(8:10),(1:2)},'label',{'fSAD','emph','PD','Norm'});
    nprm = length(prmclass);
    MF_label = {'V_e4','S_e4','B_e5','X_e5'};

    % Sub-classify PRM:
    prm_sub = zeros(size(prm));
    for i = 1:nprm
        prm_sub(ismember(prm,prmclass(i).val)) = i;
    end
    
    % Generate MF values
    fprintf('Calculating Minkowski functionals ...\n');
    p = minkowskiFun(seg,'thresh',1:nprm,'tmode','==','n',[10,10,10],'gridsp',[5,5,5],'voxsz',voxsz,'mask',seg);
    
    % Re-grid, scale, and save as .nii.gz
    scale=[10000 10000 1e5 1e5];
    for iprm = 1:nprm
        for imf = 1:4
            fn_out = [fn_base,'_',prmclass(iprm).label,'_',MF_label{imf},'.nii'];
            fprintf('Saving MF result: %s\n',fn_out);
            tprm = grid2img(p.MF(iprm,imf,:),p.ind,seg,3,1) * scale(imf);
            niftiwrite(int16(tprm),fullfile(sv_path,fn_out),info_sv,'Compressed',true);
        end
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