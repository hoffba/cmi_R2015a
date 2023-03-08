function res = CTlung_Pipeline_sub(procdir,opts)
% Pipeline processing on Great Lakes

import mlreportgen.report.*

tt = tic;

try

    fn_log = fullfile(procdir,'pipeline_log.txt');
    
    if strcmp(opts.cluster,'GL')
        opts.save_path = checkTurboPath(opts.save_path);
        opts.save_path = opts.save_path{1};
    end
    
%     R = Report('output','pdf');
    
    QC_nslice = 25;
    fn_ext = '.nii.gz';

    nowstr = datestr(datetime('now'),0);
    writeLog(fn_log,'\n\nStarting Pipeline Process : %s\n',nowstr);
    writeLog(fn_log,'GIF location: %s\n',opts.save_path);
    
    % Initialize results table
    [~,ID] = fileparts(procdir);
    regionnames = {'WholeLung','RL','LL','RUL','RML','RULplus','RLL','LUL','LULplus','LLi','LLL'};
    Nr = numel(regionnames);
    fn_res = fullfile(procdir,sprintf('%s_PipelineResults.csv',ID));
    if exist(fn_res,'file')
        res = readtable(fn_res,'Delimiter',',');
    else
        varnames = {'ID','Exp_Source','Ins_Source'};
        Nv = numel(varnames);
        res = table('Size',[Nr,Nv],'VariableTypes',repmat({'cellstr'},1,Nv),'VariableNames',varnames);
        res.ID(:) = {ID};
        
        % Read source paths from file
        tname = fullfile(procdir,sprintf('%s_SourceData.csv',ID));
        if exist(tname,'file')
            writeLog(fn_log,'Reading source location(s) from file: %s\n',tname);
            T_opts = detectImportOptions(tname,'Delimiter',',');
            T = readtable(tname,T_opts);
            res.Exp_Source(:) = T.DataPath(find(strcmp(T.Tag,'Exp'),1));
            res.Ins_Source(:) = T.DataPath(find(strcmp(T.Tag,'Ins'),1));
        else
            writeLog(fn_log,'Source files not found: %s\n',tname);
        end
    end
    res.Properties.RowNames = regionnames;
    res.ROI = regionnames';
    fld = {'ID','Exp_Source','Ins_Source'};
    for i = 1:numel(fld)
        if isnumeric(res.(fld{i}))
            res.(fld{i}) = arrayfun(@num2str,res.(fld{i}),'UniformOutput',false);
        end
    end

    % Load images
    img = struct('flag',{false,false},'mat',{[],[]},'info',{[],[]},'label',{[],[]},'QCind',{[],[]});
    fn_exp = fullfile(procdir,[ID,'.exp',fn_ext]);
    if exist(fn_exp,'file')
        writeLog(fn_log,'Reading EXP image ...\n');
        [img(1).mat,label,fov,orient,info] = cmi_load(1,[],fn_exp);
        d = size(img(1).mat);
        voxsz = fov./d;
        img(1).info = struct('label',label,'fov',fov,'orient',orient,'d',d,'voxsz',voxsz,'natinfo',info,...
            'voxvol',prod(voxsz),'name',label);
        img(1).flag = true;
    end
    fn_ins = fullfile(procdir,[ID,'.ins',fn_ext]);
    if exist(fn_ins,'file')
        writeLog(fn_log,'Reading INS image ...\n');
        [img(2).mat,label,fov,orient,info] = cmi_load(1,[],fn_ins);
        d = size(img(2).mat);
        voxsz = fov./d;
        img(2).info = struct('label',label,'fov',fov,'orient',orient,'d',d,'voxsz',voxsz,'natinfo',info,...
            'voxvol',prod(voxsz),'name',label);
        img(2).flag = true;
    end

    % Load segmentations
    fn_exp_seg = fullfile(procdir,[ID,'.exp.label',fn_ext]);
    if img(1).flag
        if exist(fn_exp_seg,'file')
            img(1).label = cmi_load(1,[],fn_exp_seg);
        else
            img(1).flag = false;
        end
    end
    fn_ins_seg = fullfile(procdir,[ID,'.ins.label',fn_ext]);
    if img(2).flag
        if exist(fn_ins_seg,'file')
            img(2).label = cmi_load(1,[],fn_ins_seg);
        else
            img(2).flag = false;
        end
    end

    % QC segmentation
    if img(1).flag
        writeLog(fn_log,'Generating EXP montage...\n');
        ind = find(std(single(img(1).mat),1,[1 2])~=0)';
        Ni = numel(ind);
        QCns = min(QC_nslice,Ni);
        if Ni == img(1).info.d(3)
            QCpad = round(img(1).info.d(3)/QC_nslice);
            img(1).QCind = round(linspace(0,img(1).info.d(3)-QCpad,QCns)+ceil(QCpad/2));
        else
            img(1).QCind = ind(round(linspace(1,Ni,QCns)));
        end
        cdata = QCmontage('seg',img(1).mat(:,:,img(1).QCind),...
                                logical(img(1).label(:,:,img(1).QCind)),...
                                img(1).info.voxsz,...
                                fullfile(procdir,sprintf('%s_Montage',img(1).info.label)),...
                                fullfile(opts.save_path,'Montage_exp.gif'));
    end
    if img(2).flag
        writeLog(fn_log,'Generating INSP montage...\n');
        ind = find(std(single(img(2).mat),1,[1 2])~=0)';
        Ni = numel(ind);
        QCns = min(QC_nslice,Ni);
        if Ni == img(2).info.d(3)
            QCpad = round(img(2).info.d(3)/QC_nslice);
            img(2).QCind = round(linspace(0,img(2).info.d(3)-QCpad,QCns)+ceil(QCpad/2));
        else
            img(2).QCind = ind(round(linspace(1,Ni,QCns)));
        end
        cdata = QCmontage('seg',img(2).mat(:,:,img(2).QCind),...
                                logical(img(2).label(:,:,img(2).QCind)),...
                                img(2).info.voxsz,...
                                fullfile(procdir,sprintf('%s_Montage',img(2).info.label)),...
                                fullfile(opts.save_path,'Montage_ins.gif'));
    end

    % Airways
    if opts.airway
            if img(2).flag
                ydir = fullfile(procdir,['yacta_',ID,'_Ins']);
                writeLog(fn_log,'YACTA directory: %s\n',ydir);
                airway_res = readYACTAairways(ydir);
                if ~isempty(airway_res)
                    
                    if isfield(airway_res,'Wall_pct__3_8_')
                        res = addTableVarVal(res,'WallPct_3_8','WholeLung',airway_res.Wall_pct__3_8_);
                    end

                    fld = {'Wall_pct',  'Wall_pct',                     'WholeLung';...
                           'Wall_pct',  'Wall_pct_RightUpperLobe',      'RUL';...
                           'Wall_pct',  'Wall_pct_RightMidLobe',        'RML';...
                           'Wall_pct',  'Wall_pct_RightUpperLobePlus',  'RULplus';...
                           'Wall_pct',  'Wall_pct_RightLowerLobe',      'RLL';...
                           'Wall_pct',  'Wall_pct_LeftUpperLobe',       'LUL';...
                           'Wall_pct',  'Wall_pct_LeftUpperLobePlus',   'LULplus';...
                           'Wall_pct',  'Wall_pct_LeftLingula',         'LLi';...
                           'Wall_pct',  'Wall_pct_LeftLowerLobe',       'LLL';
                           'Pi10',      'Pi10',         'WholeLung';...
                           'Pi10',      'Pi10_RUL',     'RUL';...
                           'Pi10',      'Pi10_RML',     'RML';...
                           'Pi10',      'Pi10_RULplus', 'RULplus';...
                           'Pi10',      'Pi10_RLL',     'RLL';...
                           'Pi10',      'Pi10_LUL'      'LUL';...
                           'Pi10',      'Pi10_LULplus'  'RULplus';...
                           'Pi10',      'Pi10_LLi',     'LLi';...
                           'Pi10',      'Pi10_LLL',     'LLL';...
                           'Pi15',      'Pi15',         'WholeLung';...
                           'BEI',       'BEI_Lung',     'WholeLung';...
                           'BEI',       'BEI_Right',    'RL';...
                           'BEI',       'BEI_Left',     'LL';...
                           'BEI',       'BEI_RUL',      'RUL';...
                           'BEI',       'BEI_RML',      'RML';...
                           'BEI',       'BEI_RULplus',  'RULplus';...
                           'BEI',       'BEI_RLL',      'RLL';...
                           'BEI',       'BEI_LUL',      'LUL';...
                           'BEI',       'BEI_LULplus',  'LULplus';...
                           'BEI',       'BEI_LLi',      'LLi';...
                           'BEI',       'BEI_LLL',      'LLL';...
                           'BEI_gen',   'x1_Level_Lung',    'WholeLung';...
                           'BEI_gen',   'x1_Level_Right',   'RL';...
                           'BEI_gen',   'x1_Level_Left',    'LL';...
                           'BEI_gen',   'x1_Level_RUL',     'RUL';...
                           'BEI_gen',   'x1_Level_RML',     'RML';...
                           'BEI_gen',   'x1_Level_RULplus', 'RULplus';...
                           'BEI_gen',   'x1_Level_RLL',     'RLL';...
                           'BEI_gen',   'x1_Level_LUL',     'LUL';...
                           'BEI_gen',   'x1_Level_LLi',     'LLi';...
                           'BEI_gen',   'x1_Level_LULplus', 'LULplus';...
                           'BEI_gen',   'x1_Level_LLL',     'LLL'};
                    for ifld = 1:size(fld,1)
                        if isfield(airway_res,fld{ifld,2})
                            res = addTableVarVal(res,fld{ifld,1},fld{ifld,3},airway_res.(fld{ifld,2}));
                        end
                    end

                    genstr = {'WT'};
                    segstr = { 'seg',   4   ;...
                               'subseg',5:7 };
                    lobestr =  {'Right','Left','RUL','RML','RUL+','RLL','LUL','LUL+','LLi','LLL'};
                    lobestr2 = {'RL','LL','RUL','RML','RULplus','RLL','LUL','LULplus','LLi','LLL'};
                    for i = 1:numel(genstr)
                        for j = 1:2
                            for k = 1:numel(lobestr)
                                if isfield(airway_res,genstr{i}) && isfield(airway_res.(genstr{i}),lobestr{k})
                                    segind = segstr{j,2};
                                    vals = airway_res.(genstr{i}).(lobestr{k});
                                    vals = vals(segind(segind<=numel(vals)));
                                    res = addTableVarVal(res,[genstr{i},'_',segstr{j,1}],lobestr2{k},mean(vals(vals>0)));
                                end
                            end
                        end
                    end

                end
            end
    end

    % ScatterNet for AT on Exp CT scan
    atMap = [];
    try
        if opts.scatnetAT && img(1).flag
            fn_scatnet = fullfile(procdir,sprintf('%s.%s%s',res.ID{1},'scatnet_AT',fn_ext));
            writeLog(fn_log,'Air trapping map ... ');
            if exist(fn_scatnet,'file')
                writeLog(fn_log,'from file\n');
                atMap = cmi_load(1,img(1).info.d(1:3),fn_scatnet);
            else
                writeLog(fn_log,'generating with ScatNet\n');
                atMap = ScatNet(img(1).mat,logical(img(1).label),0);
                cmi_save(0,atMap,'ScatNet',img(1).info.fov,img(1).info.orient,fn_scatnet);
            end
        end
    catch err
        writeLog(fn_log,'ScatNet FAILED:\n%s\n',getReport(err));
    end
    
    % Vessel analysis
    try
        if opts.vessel && img(2).flag
            writeLog(fn_log,'Vessel analysis ...\n');
            fn_re_ins = fullfile(procdir,sprintf('re_%s.ct.nii.gz',res.ID{1}));
            fn_re_seg = fullfile(procdir,sprintf('re_%s.lobe_segmentation.nii.gz',res.ID{1}));
            if exist(fn_re_ins,'file') && exist(fn_re_seg,'file')
                tinfo = niftiinfo(fn_re_ins);
                tins = niftiread(tinfo);
                tseg = niftiread(fn_re_seg);
            else
                tinfo = init_niftiinfo(res.ID{1},img(2).info.voxsz,class(img(2).mat),img(2).info.d);
                tins = img(2).mat;
                tseg = img(2).label;
            end
            T = pipeline_vesselseg( tins, tseg, tinfo, procdir, opts, fn_log);
            res = addTableVarVal(res,T);
        end
    catch err
        writeLog(fn_log,'Vessel Analysis FAILED:\n%s\n',getReport(err));
    end

    % Quantify unregistered CT scans
    writeLog(fn_log,'Quantifying unregistered statistics\n');
    if opts.unreg && img(1).flag
        T = CTlung_Unreg('exp',img(1).mat,img(1).info.voxvol,img(1).label,atMap);
        clear atMap;
        res = addTableVarVal(res,T);
    end
    if opts.unreg && img(2).flag
        T = CTlung_Unreg('ins',img(2).mat,img(2).info.voxvol,img(2).label);
        res = addTableVarVal(res,T);
    end

    % Register I2E
    if img(1).flag && img(2).flag && (opts.reg || opts.prm || opts.tprm || opts.jac)
        ins_reg = [];
        fn_reg = fullfile(procdir,sprintf('%s.%s%s',res.ID{1},'ins.reg',fn_ext));
        elxdir = fullfile(procdir,sprintf('elastix_%s',ID));
        if exist(fn_reg,'file')
            writeLog(fn_log,'Loading registered INS from file...\n');
            ins_reg = readNIFTI(fn_reg);
        elseif opts.reg
            writeLog(fn_log,'Performing registration...\nelxdir : %s\n',elxdir)
            pipeline_reg( img , elxdir , res.ID{1} , opts );
            
            % Move registered file out to procdir
            [~,outfn] = fileparts(img(2).info.name);
            fn_move = {'spatialJacobian',       '.jac';...
                       'FullSpatialJacobian',   '.fulljac';...
                       'DeformationFields',     '.def';...
                       [outfn,'_R'],            '.ins.reg'};
            nf = size(fn_move,1);
            for i = 1:nf
                fn = dir(fullfile(elxdir,[fn_move{i,1},'*']));
                if ~isempty(fn)
                    writeLog(fn_log,'Re-saving: %s\n',fn_move{i,1});
                    fn_elx = fullfile(elxdir,fn(1).name);
                    timg = cmi_load(1,[],fn_elx);
                    cmi_save(0,timg,fn_move{i,1},img(1).info.fov,img(1).info.orient,...
                        fullfile(procdir,sprintf('%s%s%s',res.ID{1},fn_move{i,2},fn_ext)));
                    % Delete image from elxreg_ folder:
                    delete(fn_elx);
                    if endsWith(fn_elx,'.mhd')
                        delete([fn_elx(1:end-3),'raw']);
                    end
                    if i==nf
                        ins_reg = timg;
                        clear timg;
                    end
                end
            end
        end

        % QC registration
        writeLog(fn_log,'Saving Registration Montage ...\n');
        QCmontage('reg',ins_reg(:,:,img(1).QCind),...
                        logical(img(1).label(:,:,img(1).QCind)),...
                        img(1).info.voxsz,...
                        fullfile(procdir,sprintf('%s_Reg_Montage',res.ID{1})),...
                        fullfile(opts.save_path,'Montage_reg.gif'));
        img(2) = [];

        % Jacobian analysis
        fn_jac = fullfile(procdir,sprintf('%s.jac%s',res.ID{1},fn_ext));
        if opts.jac
            jac = [];
            if exist(fn_jac,'file')
                writeLog(fn_log,'Analyzing Jacobian map...\n');
                jac = cmi_load(1,[],fn_jac);
            else
                writeLog(fn_log,'Jacobian file not found.\n');
                if isfolder(elxdir)
                    writeLog(fn_log,'Calculating Jacobian from Transform Parameters.\n');
                    
                    fn = dir(fullfile(elxdir,'TransformParameters.?.txt'));
                    fn = fullfile(elxdir,fn(end).name);
                    
                    % Fix transform parameter filenames in the transform chain
                    fixTransformParameter(fn);
                    
                    elxObj = ElxClass;
                    str = elxObj.sysCmd(elxdir,'tp',fn,'jac',true,'wait',true);
                    system(str);
                    clear elxObj;
                    
                    % Re-save Jacobian map to subject folder:
                    fn_jac_elx = fullfile(elxdir,'spatialJacobian.nii');
                    if exist(fn_jac_elx,'file')
                        jac = cmi_load(1,[],fn_jac_elx);
                        cmi_save(0,jac,'spatialJacobian',img(1).info.fov,img(1).info.orient,fn_jac);
                        delete(fn_jac_elx);
                    else
                        writeLog(fn_log,'Jacobian processing failed.\n');
                    end
                end
            end
            if ~isempty(jac)
                T = lobeLoop(img.label,@(mask,A,str)tabulateStats(mask,A,str),jac,'Jac');
                res = addTableVarVal(res,T);
            end
        end
        
        % ScatNet for Emph
        SNemph = [];
        try
            if opts.scatnetEmph && img(1).flag
                fn_SNemph = fullfile(procdir,sprintf('%s.%s%s',res.ID{1},'scatnet_Emph',fn_ext));
                writeLog(fn_log,'Air trapping map ... ');
                if exist(fn_SNemph,'file')
                    writeLog(fn_log,'from file\n');
                    SNemph = cmi_load(1,img(1).info.d(1:3),fn_SNemph);
                else
                    writeLog(fn_log,'generating with ScatNet\n');
                    SNemph = ScatterNetEMPH(img(1).mat,logical(img(1).label),0);
                    cmi_save(0,SNemph,'ScatNet',img(1).info.fov,img(1).info.orient,fn_SNemph);
                end
            end
        catch err
            writeLog(fn_log,'ScatNet FAILED:\n%s\n',getReport(err));
        end
        
        % Blood density change map:
        if opts.dBlood
            dBlood = [];
            fn_dBlood = fullfile(procdir,sprintf('%s.%s%s',res.ID{1},'dblood',fn_ext));
            if exist(fn_dBlood,'file')
                writeLog(fn_log,'Loading dBlood from file...\n');
                dBlood = readNIFTI(fn_dBlood);
            elseif exist('ins_reg','var') && exist('jac','var')
                writeLog(fn_log,'Calculating dBlood...\n');
                dBlood = pipeline_blood_density(img(1).mat,ins_reg,jac,img(1).label);
                cmi_save(0,dBlood,{'dBlood'},img(1).info.fov,img(1).info.orient,fn_dBlood);
            end
            if ~isempty(dBlood)
                T = lobeLoop(img.label,@(mask,A,str)tabulateStats(mask,A,str),dBlood,'dBlood');
                res = addTableVarVal(res,T);
            end
        end
        clear jac
        
        % PRM calculation
        prm10 = [];
        if opts.prm
            fn_PRM = fullfile(procdir,sprintf('%s.%s%s',res.ID{1},'prm',fn_ext));
            if exist(fn_PRM,'file')
                writeLog(fn_log,'Loading PRM from file...\n');
                prm10 = int8(readNIFTI(fn_PRM));
            else
                writeLog(fn_log,'Calculating PRM...\n');
                [prm10,~] = pipeline_PRM(img(1).mat,img(1).info,logical(img(1).label),ins_reg,...
                    fullfile(procdir,sprintf('%s_PRM_Scatter',res.ID{1})),...
                    fullfile(opts.save_path,'Montage_PRM_Scatter.gif'));

                % Save PRM
                writeLog(fn_log,'Saving PRM as NIFTI ... ');
                stat = cmi_save(0,prm10,{'PRM'},img.info.fov,img.info.orient,fn_PRM);
                if stat
                    writeLog(fn_log,'  PRM saved\n');
                else
                    writeLog(fn_log,'  Could not save PRM to file.\n');
                end
            end
        end
        clear ins_reg;

        if ~isempty(prm10)
            % Tabulate 10-color PRM results
            writeLog(fn_log,'Tabulating 10-color PRM results...\n');
            T = lobeLoop(img.label,@(mask,prm,flag)tabulatePRM(mask,prm,flag),prm10,1);
            res = addTableVarVal(res,T);

            % map full PRM values (1:10) to (norm,fsad,emph,pd,ns)
            prmlabel = {'Norm', 'fSAD', 'Emph', 'PD',       'NS';...
                        [1,2],  3,      [4,5],  [8,9,10],   6    };
            prm5 = int8(zeros(size(prm10)));
            for i = 1:size(prmlabel,2)
                prm5(ismember(prm10,prmlabel{2,i})) = i;
            end
            clear prm10

            % QC PRM
            writeLog(fn_log,'Generating PRM Montage ...\n');
            QCmontage('prm',img.mat(:,:,img(1).QCind),...
                            double(prm5(:,:,img(1).QCind)),...
                            img.info.voxsz,...
                            fullfile(procdir,sprintf('%s_PRM_Montage',res.ID{1})),...
                            fullfile(opts.save_path,'Montage_PRM.gif'));

            % Tabulate 5-color PRM results
            writeLog(fn_log,'Tabulating 5-color PRM results...\n');
            T = lobeLoop(img.label,@(mask,prm,flag)tabulatePRM(mask,prm,flag),prm5,0);
            res = addTableVarVal(res,T);

            % Calculate tPRM
            if opts.tprm
                prmlabel = ["norm","fsad","emph","pd"];
                mflabel = ["v","s","b","x"];
                fn_tprm = fullfile(procdir,...
                    string(res.ID{1})+".tprm."+prmlabel+"."+mflabel'+string(fn_ext));
                if all(cellfun(@(x)exist(x,'file'),fn_tprm))
                    writeLog(fn_log,'Loading tPRM from files ...\n');
                    clear prm5;
                    for iprm = 1:numel(prmlabel)
                        for imf = 1:numel(mflabel)
                            writeLog(fn_log,'   %s - %s\n',prmlabel(iprm),mflabel(imf))
                            tprm = readNIFTI(fn_tprm(imf,iprm));
                            str = prmlabel(iprm)+'_'+upper(mflabel(imf));
                            T = lobeLoop(img.label,@(mask,tprm,str)tabulateTPRM(mask,tprm,str),tprm,str);
                            res = addTableVarVal(res,T);
                        end
                    end
                elseif opts.tprm
                    t = tic;
                    writeLog(fn_log,'Generating tPRM ...\n');

                    % Calculate MF values
                    p = minkowskiFun(prm5,'thresh',1:4,...
                        'tmode','==',...
                        'n',10*ones(1,3),...
                        'gridsp',5,...
                        'voxsz',img.info.voxsz,...
                        'mask',logical(img.label),...
                        'prog',0);
                    clear prm5;

                    % Interpolate to maps
                    for ithresh = 1:size(p.MF,1)
                        for imf = 1:size(p.MF,2)
                            writeLog(fn_log,'   %s - %s\n',prmlabel(ithresh),mflabel(imf));

                            % Interpolate to image space
                            tstr = prmlabel(ithresh) + '.' + mflabel(imf);
                            writeLog(fn_log,'       Interpolating\n');
                            tprm = grid2img(p.MF(ithresh,imf,:),p.ind,p.mask,3,1);

                            % Save tPRM image
                            writeLog(fn_log,'       Saving NIFTI\n');
                            cmi_save(0,single(tprm),{char(tstr)},img.info.fov,img.info.orient,char(fn_tprm(imf,ithresh)));

                            % Tabulate statistics
                            writeLog(fn_log,'       Tabulating means\n');
                            str = prmlabel(ithresh)+'_'+upper(mflabel(imf));
                            T = lobeLoop(img.label,@(mask,tprm,str)tabulateTPRM(mask,tprm,str),tprm,str);
                            res = addTableVarVal(res,T);
                        end
                    end
                    writeLog(fn_log,'... tPRM complete (%s)\n',datestr(duration(0,0,toc(t)),'HH:MM:SS'));

                end
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
    writeLog(fn_log,'ERROR in CTlung_Pipeline_sub:\n%s',getReport(err,'extended','hyperlinks','off'));
    assignin('base','PipelineError',err);
end

% Need to remove rownames for future concatenation
res.Properties.RowNames = {};
    
writeLog(fn_log,'Pipeline total time = %s\n',datestr(duration(0,0,toc(tt)),'HH:MM:SS'))

function T = tabulateStats(mask,A,str)
    vname = strcat(str,{'_mean','_var'});
    nv = numel(vname);
    T = table('Size',[1,nv],'VariableTypes',repmat({'double'},1,nv),'VariableNames',vname);
    T.(vname{1}) = mean(A(mask));
    T.(vname{2}) = var(A(mask));

function T = tabulatePRM(mask,prm,flag)
    if flag % 10-color
        vals = 1:10;
        tag = cellfun(@num2str,num2cell(vals),'UniformOutput',false);
    else % 5-color
        vals = 1:5;
        tag = {'Norm', 'fSAD', 'Emph', 'PD', 'NS'};
    end
    nv = numel(vals);
    vname = strcat('PRM_',tag);
    T = table('Size',[1,nv],'VariableTypes',repmat({'double'},1,nv),...
        'VariableNames',vname);
    np = nnz(mask);
    for i = 1:numel(vals)
        T.(vname{i}) = nnz(prm(mask)==i)/np*100;
    end
    
function T = tabulateTPRM(mask,tprm,str)
    vname = sprintf('tPRM_%s',str);
    T = table('Size',[1,1],'VariableTypes',{'double'},'VariableNames',{vname});
    T.(vname) = mean(tprm(mask));
    
function img = medfilt2_3(img)
    for i = 1:size(img,3)
        img(:,:,i) = medfilt2(img(:,:,i));
    end
   