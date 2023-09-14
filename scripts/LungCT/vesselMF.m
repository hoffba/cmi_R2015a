function res = vesselMF(varargin)
% Performs Minkowski functional analysis on vessel maps
% Inputs (GL batch): vesselMF(procdir,method,username)
%   procdir     = directory containing all data
%   method      = <'GL','batch','parbatch'>
%   username    = UMich uniquename
% Inputs (single case): vesselMF(true,ID,procdir,fn_results)
%   ID          = string identification for case
%   procdir     = directory of data
%   fn_results  = filename for results .xlsx

    res = [];
    if nargin==3 && islogical(varargin{1})
        % Run single process (e.g. on Great Lakes, or for debugging)
        res = vesselMF_sub(varargin{2:3});
    else
        opts = parseInputs(varargin);

        fn_results = fullfile(opts.procdir,...
            sprintf('Minkowski_output_%s.xlsx',char(datetime('now','Format','yyyyMMddHHmmss'))));
        
        % Find data to process
        % ** may need to adjust input parameters based on file structure
        % Files all saved in same folder
        fn = dir(fullfile(opts.procdir,'*.ct.nii.gz'));
        ID = cellfun(@(x)extractBefore(x,'.'),{fn.name}','UniformOutput',false);
        ID(startsWith(ID,'re_')) = [];
        
        switch opts.method
            case 'debug'
                batchloop(ID,opts.procdir,fn_results);
            case 'batch'
                batch(@batchloop,0,{ID,opts.procdir,fn_results});
            case 'parbatch'
                batch(@parbatchloop,0,{ID,opts.procdir,fn_results});
            case 'GL'
                GL_run(opts.username, 'vesselMF',...
                    {true,   ID,     opts.procdir},...
                    [false,  false,  true],... % flag for path names
                    [true,   false,  true],...  % flag for static inputs
                    'ProcessMemory',opts.mem,'ProcessTime',opts.time,...
                    'Partition',opts.part,'Nodes',opts.nnodes,...
                    'TimeStamp',char(datetime('now','Format','yyyyMMddHHmmss')),...
                    'save_path',opts.procdir);
        end
    end

function opts = parseInputs(var_in)
    p = inputParser;
    addParameter(p,'procdir',pwd,@isfolder);
    addParameter(p,'method','parbatch',@(x)ismember(x,{'debug','batch','parbatch','GL'}));
    addParameter(p,'username','',@ischar);
    addParameter(p,'mem',24,@isnumeric);
    addParameter(p,'time',200,@isnumeric);
    addParameter(p,'part','standard',@(x)ismember(x,{'standard','largemem'}));
    addOptional(p,'nnodes',1,@isnumeric);
    parse(p,var_in{:});
    opts = p.Results;

function batchloop(ID,procdir)
    N = numel(ID);
    for i = 1:N
        fprintf('\nStarting process #%d of %d: %s\n',i,N,ID{i});
        t = tic;
        vesselMF_sub(ID{i},procdir);
        
        fprintf('... Finished (%s)\n',duration(0,0,toc(t)));
    end
    
function parbatchloop(ID,procdir,fn_results)
    T = [];
    N = numel(ID);
    f(1:N) = parallel.FevalFuture;
    for i = 1:N
        fprintf('Starting process #%d of %d: %s\n',i,N,ID{i});
    	f(i) = parfeval(@(x,y)vesselMF_sub(x,y),1,ID{i},procdir);
    end
    % Flag processes as they complete
    for i = 1:N
        try 
            [idx,val] = fetchNext(f);
            fprintf('Finished process #%d of %d: %s\n', idx, N, ID{i});
            if istable(val) && ~isempty(val)
                val.Properties.RowNames = {};
                T = [T;val];
            end
        catch err
            fprintf('%s\n',getReport(err,'extended','hyperlinks','off'));
        end
    end
    if ~isempty(T)
        writetable(T,fn_results);
    end

function res = vesselMF_sub(ID,procdir)

    res = [];

    fn_results = fullfile(procdir,sprintf('%s_Results.csv',ID));
    if isfile(fn_results)
        fprintf('Reading results from saved table: %s\n',fn_results)
        res = readtable(fn_results);
        return;
    end

    fprintf('Performing vessel MF analysis: %s\n',ID);

    fn_ct = fullfile(procdir,[ID,'.ct.nii.gz']);
    fn_lung = dir(fullfile(procdir,[ID,'.*.lobe_segmentation.nii.gz']));
    if isempty(fn_lung)
        fprintf('   Segmentation not found\n')
        return;
    end
    fn_lung = fullfile(procdir,fn_lung.name);
    if ~(isfile(fn_ct) && isfile(fn_lung))
        return;
    end

    fn_re_lung = fullfile(procdir,['re_',ID,'.lobe_segmentation.nii.gz']);
    fn_vessels = fullfile(procdir,[ID,'_binVessel.nii.gz']);

    fn_vess_results = fullfile(procdir,[ID,'_vesselMetrics.csv']);
    if isfile(fn_vessels) && isfile(fn_vess_results) && isfile(fn_re_lung)
        % Load vessel results
        res = readtable(fn_vess_results);
    else
        % Generate vessel segmentation
        info = niftiinfo(fn_ct);
        img = single(niftiread(info));
        seg_lung = niftiread(fn_lung);
        opts = struct('frangi',true,'curvi',true);
        info.Description = ID;
        res = pipeline_vesselseg( img, seg_lung, info, procdir, opts, []);
    end

    % Add ID to table
    res = addvars(res,repmat({ID},size(res,1),1),'Before',1,'NewVariableNames','ID');

    % Load resampled lung segmentation
    info = niftiinfo(fn_re_lung);
    seg_lung = niftiread(info);
        
    % Minkowski Functional map file names
    fnV = fullfile(procdir,[ID,'.vessel.v.nii.gz']);
    fnS = fullfile(procdir,[ID,'.vessel.s.nii.gz']);
    fnB = fullfile(procdir,[ID,'.vessel.b.nii.gz']);
    fnX = fullfile(procdir,[ID,'.vessel.x.nii.gz']);
    
    if isfile(fnV) && isfile(fnS) && isfile(fnB) && isfile(fnX)
        % Load maps from file
        fprintf('Reading MF maps from file...\n')
        MF(:,:,:,1) = niftiread(fnV);
        MF(:,:,:,2) = niftiread(fnS);
        MF(:,:,:,3) = niftiread(fnB);
        MF(:,:,:,4) = niftiread(fnX);
    else
        % Load binarized vessel map
        seg_vess = niftiread(fn_vessels)==1;

        % Calculate Minkowski functionals
        fprintf('Calculating vessel MF maps...\n');
        p = minkowskiFun(seg_vess,'n',[10,10,10],'gridsp',[5,5,5],'voxsz',ones(1,3),'mask',seg_lung>0);
        MF = single(grid2img(p.MF(1,:,:),p.ind,p.mask,3,1));

        % Save maps to Nifti
        fprintf('Saving MF maps...\n');
        info.Datatype = 'single';
        niftiwrite(MF(:,:,:,1),fnV(1:end-3),info,'Compressed',true);
        niftiwrite(MF(:,:,:,2),fnS(1:end-3),info,'Compressed',true);
        niftiwrite(MF(:,:,:,3),fnB(1:end-3),info,'Compressed',true);
        niftiwrite(MF(:,:,:,4),fnX(1:end-3),info,'Compressed',true);
    end
    
    % Print slices in JPEG
%     saveJPG(seg_lung,seg_vess,MF);

    % Add MF metrics to results table
    res =  addTableVarVal(res,CTlung_LobeStats(seg_lung,'Vessel_V','mean',MF(:,:,:,1)));
    res =  addTableVarVal(res,CTlung_LobeStats(seg_lung,'Vessel_S','mean',MF(:,:,:,2)));
    res =  addTableVarVal(res,CTlung_LobeStats(seg_lung,'Vessel_B','mean',MF(:,:,:,3)));
    res =  addTableVarVal(res,CTlung_LobeStats(seg_lung,'Vessel_X','mean',MF(:,:,:,4)));

    % Write results to Excel file
    writetable(res,fn_results);

function saveJPG(seg_lung,seg_vess,MF)
    ns = size(seg_vess,3);
    ind = any(seg_lung,[1,2]);
    fn_jpg = fullfile(procdir,arrayfun(@(x)sprintf([ID,'_%04d.jpg'],x),1:ns,'UniformOutput',false)');
    if ~all(isfile(fn_jpg(ind)))
        % Set up figure
        hf = figure('visible','off');
        ha1 = subplot(2,3,1,'Parent',hf);
            him1 = imshow(seg_vess(:,:,1),[0,1],'Parent',ha1);
            title(ha1,'vessel\_seg', 'FontSize', 15);
            colorbar(ha1);
        ha2 = subplot(2,3,2,'Parent',hf);
            him2 = imshow(seg_lung(:,:,1),[0,1],'Parent',ha2);
            title(ha2,'lung\_seg', 'FontSize', 15);
            colorbar(ha2);
        ha3 = subplot(2,3,3,'Parent',hf);
            him3 = imshow(MF(:,:,1,1),[0,.6],'Parent',ha3);
            title(ha3,'V', 'FontSize', 15);
            colorbar(ha3);
        ha4 = subplot(2,3,4,'Parent',hf);
            him4 = imshow(MF(:,:,1,1),[0,.3],'Parent',ha4);
            title(ha4,'S', 'FontSize', 15);
            colorbar(ha4);
        ha5 = subplot(2,3,5,'Parent',hf);
            him5 = imshow(MF(:,:,1,1),[0,.015],'Parent',ha5);
            title(ha5,'B', 'FontSize', 15);
            colorbar(ha5);
        ha6 = subplot(2,3,6,'Parent',hf);
            him6 = imshow(MF(:,:,1,1),[-.005,.005],'Parent',ha6);
            title(ha6,'Chi', 'FontSize', 15);
            colorbar(ha6);

        % Loop over slices
        fprintf('Saving slices to JPEG...\n');
        for i = 1:size(seg_vess,3)
            if ind(i)
                try
                    him1.CData = seg_vess(:,:,i);
                    him2.CData = seg_lung(:,:,i);
                    him3.CData = MF(:,:,i,1);
                    him4.CData = MF(:,:,i,2);
                    him5.CData = MF(:,:,i,3);
                    him6.CData = MF(:,:,i,4);
                    saveas(hf, fn_jpg{i});
                    pause(0.001);
                catch err
                    fprintf('ERROR: slice #%d\n%s\n',i,getReport(err,'extended','hyperlinks','off'));
                end
            end
        end
        delete(hf);
    end
