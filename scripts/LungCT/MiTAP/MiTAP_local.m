function res = MiTAP_local(procdir,opts,varargin)
% Pipeline processing on Gre

tt = tic;
res = table;
try

    % Determine inputs
    [~,ID] = fileparts(procdir);
    opts.exp_src = '';
    opts.ins_src = '';
    if nargin==4
        opts.exp_src = varargin{1};
        opts.ins_src = varargin{2};

        % Save data sources to table in processing directory for future reference
        T = table('Size',[2,2],'VariableTypes',{'cellstr','cellstr'},'VariableNames',{'Tag','DataPath'});
        T.Tag(1) = {'Exp'};
        T.DataPath(1) = {exp_src};
        T.Tag(2) = {'Ins'};
        T.DataPath(2) = {ins_src};
        writetable(T,fullfile(procdir,sprintf('%s_SourceData.csv',ID)));
    end

    % Initialize case-specific options
    opts.procdir = procdir;
    opts.ID = ID;
    fn_ext = '.nii.gz';
    % ~~ Log file ~~
    opts.fn_log = fullfile(procdir,'pipeline_log.txt');
    % ~~ Results CSV file
    fn_res = fullfile(procdir,[ID,'_PipelineResults.csv']);
    % tPRM labels
    prmlabel = ["norm","fsad","emph","pd","ns"];
    mflabel = ["v","s","b","x"];
    % ~~ Image files ~~
    opts.fn.exp =               fullfile(procdir,[ID,'.exp',fn_ext]);
    opts.fn.ins =               fullfile(procdir,[ID,'.ins',fn_ext]);
    opts.fn.exp_seg =           fullfile(procdir,[ID,'.exp.label',fn_ext]);
    opts.fn.ins_seg =           fullfile(procdir,[ID,'.ins.label',fn_ext]);
    opts.fn.scatnetAT =         fullfile(procdir,[ID,'.scatnetAT',fn_ext]);
    opts.fn.scatnetAT_PEDS =    fullfile(procdir,[ID,'.scatnetAT_PEDS',fn_ext]);
    opts.fn.scatnetEmph =       fullfile(procdir,[ID,'.scatnetEmph',fn_ext]);
    opts.fn.reg =               fullfile(procdir,[ID,'.ins.reg',fn_ext]);
    opts.fn.jac =               fullfile(procdir,[ID,'.jac',fn_ext]);
    opts.fn.fulljac =           fullfile(procdir,[ID,'.fulljac',fn_ext]);
    opts.fn.def =               fullfile(procdir,[ID,'.def',fn_ext]);
    opts.fn.scatnetEmphReg =    fullfile(procdir,[ID,'.scatnetEmphInsR',fn_ext]);
    opts.fn.dBlood =            fullfile(procdir,[ID,'.dblood',fn_ext]);
    opts.fn.prm =               fullfile(procdir,[ID,'.prm',fn_ext]);
    opts.fn.tprm = fullfile(procdir, string(ID) + ".tprm." + prmlabel(1:4) + "." + mflabel' + string(fn_ext));
    % ~~ QC figure files ~~
    opts.QCnslice = 25;
    opts.fn.expMontage =        fullfile(procdir,[ID,'_Exp_Montage.tif']);
    opts.fn.insMontage =        fullfile(procdir,[ID,'_Ins_Montage.tif']);
    opts.fn.regMontage =        fullfile(procdir,[ID,'_InsReg_Montage.tif']);
    opts.fn.prmMontage =        fullfile(procdir,[ID,'_PRM_Montage.tif']);
    opts.fn.prmScatter =        fullfile(procdir,[ID,'_PRM_Scatter.tif']);
    
    % Make sure GL has the correct path
    if strcmp(opts.proc.cluster,'GL')
        opts.save_path = checkTurboPath(opts.save_path);
        opts.save_path = opts.save_path{1};
    end

    % Start of Pipeline processing
    writeLog(opts.fn_log,'\n\nStarting Local Pipeline Process : %s\n',datetime('now'));
    writeLog(opts.fn_log,'Report location: %s\n',opts.save_path);
    
    % Initialize results table
    regionnames = {'WholeLung','RL','LL','RUL','RML','RULplus','RLL','LUL','LULplus','LLi','LLL'};
    Nr = numel(regionnames);
    if exist(fn_res,'file')
        res = readtable(fn_res,'Delimiter',',');
        if isnumeric(res.ID)
            res.ID = arrayfun(@num2str,res.ID,'UniformOutput',false);
        end
    else
        varnames = {'ID','ROI'};
        Nv = numel(varnames);
        res = table('Size',[Nr,Nv],'VariableTypes',repmat({'cellstr'},1,Nv),'VariableNames',varnames);
        res.ID(:) = {ID};
        res.ROI = regionnames';
    end
    res.Properties.RowNames = regionnames;

    % Initialize image structure
    img = MiTAP_LoadImages(opts);

    % Loop over modules that need to be run locally
    for i = 1:numel(opts.mods)
        if opts.mods(i).local
            try
                [img,res] = feval(opts.modules{i},img,opts,res);
            catch err
                writeLog(opts.fn_log,'Module %s FAILED:\n%s\n',opts.modules{i},getReport(err))
            end
        end
    end

    % Save Results Table:
    if istable(res)
        writetable(res,fn_res);
    else
        fn_res = regexprep(fn_res,'.csv','.mat');
        save(fn_res,'res');
    end

catch err
    writeLog(opts.fn_log,'ERROR in MiTAP_local:\n%s',getReport(err,'extended','hyperlinks','off'));
    assignin('base','PipelineError',err);
end

% Need to remove rownames for future concatenation
res.Properties.RowNames = {};
    
writeLog(opts.fn_log,'Pipeline total time = %s\n',datetime([0,0,0,0,0,toc(tt)],'Format','HH:mm:ss'))
