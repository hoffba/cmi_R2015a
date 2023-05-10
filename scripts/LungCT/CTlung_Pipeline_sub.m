function res = CTlung_Pipeline_sub(procdir,opts)
% Pipeline processing on Great Lakes

import mlreportgen.report.*;
import mlreportgen.dom.*;

tt = tic;

try
    % Determine file names
    [~,ID] = fileparts(procdir);
    fn_ext = '.nii.gz';
    fn_log = fullfile(procdir,'pipeline_log.txt');
    fn_res = fullfile(procdir,[ID,'_PipelineResults.csv']);
    % tPRM labels
    prmlabel = ["norm","fsad","emph","pd","ns"];
    mflabel = ["v","s","b","x"];
    % ~~ Image files ~~
    opts.fn.exp = fullfile(procdir,[ID,'.exp',fn_ext]);
    opts.fn.ins = fullfile(procdir,[ID,'.ins',fn_ext]);
    opts.fn.exp_seg = fullfile(procdir,[ID,'.exp.label',fn_ext]);
    opts.fn.ins_seg = fullfile(procdir,[ID,'.ins.label',fn_ext]);
    opts.fn.scatnetAT = fullfile(procdir,[ID,'.scatnetAT',fn_ext]);
    opts.fn.scatnetEmph = fullfile(procdir,[ID,'.scatnetEmph',fn_ext]);
    opts.fn.reg = fullfile(procdir,[ID,'.ins.reg',fn_ext]);
    opts.fn.jac = fullfile(procdir,[ID,'.jac',fn_ext]);
    opts.fn.scatnetEmphReg = fullfile(procdir,[ID,'.scatnetEmphInsR',fn_ext]);
    opts.fn.dBlood = fullfile(procdir,[ID,'.dblood',fn_ext]);
    opts.fn.prm = fullfile(procdir,[ID,'.prm',fn_ext]);
    opts.fn.tprm = fullfile(procdir, string(ID) + ".tprm." + prmlabel(1:4) + "." + mflabel' + string(fn_ext));
    % ~~ QC figure files ~~
    opts.fn.expMontage = fullfile(procdir,[ID,'_Exp_Montage.tif']);
    opts.fn.insMontage = fullfile(procdir,[ID,'_Ins_Montage.tif']);
    opts.fn.regMontage = fullfile(procdir,[ID,'_InsReg_Montage.tif']);
    opts.fn.prmMontage = fullfile(procdir,[ID,'_PRM_Montage.tif']);
    opts.fn.prmScatter = fullfile(procdir,[ID,'_PRM_Scatter.tif']);
    
    % Make sure GL has the correct path
    if strcmp(opts.cluster,'GL')
        opts.save_path = checkTurboPath(opts.save_path);
        opts.save_path = opts.save_path{1};
    end
    
    % Number of slices to show in montages
    QC_nslice = 25;

    writeLog(fn_log,'\n\nStarting Pipeline Process : %s\n',datetime('now'));
    writeLog(fn_log,'GIF location: %s\n',opts.save_path);
    
    % Initialize results table
    regionnames = {'WholeLung','RL','LL','RUL','RML','RULplus','RLL','LUL','LULplus','LLi','LLL'};
    Nr = numel(regionnames);
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
    if exist(opts.fn.exp,'file')
        writeLog(fn_log,'Reading EXP image ...\n');
        [img(1).mat,label,fov,orient,info] = cmi_load(1,[],opts.fn.exp);
        d = size(img(1).mat);
        voxsz = fov./d;
        img(1).info = struct('label',label,'fov',fov,'orient',orient,'d',d,'voxsz',voxsz,'natinfo',info,...
            'voxvol',prod(voxsz),'name',label);
        img(1).flag = true;
    end
    if exist(opts.fn.ins,'file')
        writeLog(fn_log,'Reading INS image ...\n');
        [img(2).mat,label,fov,orient,info] = cmi_load(1,[],opts.fn.ins);
        d = size(img(2).mat);
        voxsz = fov./d;
        img(2).info = struct('label',label,'fov',fov,'orient',orient,'d',d,'voxsz',voxsz,'natinfo',info,...
            'voxvol',prod(voxsz),'name',label);
        img(2).flag = true;
    end

    % Load segmentations
    if img(1).flag
        if exist(opts.fn.exp_seg,'file')
            img(1).label = cmi_load(1,[],opts.fn.exp_seg);
        else
            img(1).flag = false;
        end
    end
    if img(2).flag
        if exist(opts.fn.ins_seg,'file')
            img(2).label = cmi_load(1,[],opts.fn.ins_seg);
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
        QCmontage('seg',img(1).mat(:,:,img(1).QCind),...
                  logical(img(1).label(:,:,img(1).QCind)),...
                  img(1).info.voxsz,...
                  opts.fn.expMontage,...
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
        QCmontage('seg',img(2).mat(:,:,img(2).QCind),...
                  logical(img(2).label(:,:,img(2).QCind)),...
                  img(2).info.voxsz,...
                  opts.fn.insMontage,...
                  fullfile(opts.save_path,'Montage_ins.gif'));
    end

    % Airways
    if opts.airway
        if img(2).flag
            ydir = fullfile(procdir,['yacta_',ID,'_Ins']);
            writeLog(fn_log,'YACTA directory: %s\n',ydir);
            airway_res = readYACTAairways(ydir);
            if ~isempty(airway_res)
                % Initialize table
                T = table('Size',[Nr,1],'VariableTypes',{'cellstr'},'VariableNames',{'ROI'});
                T.ROI = regionnames';
                T.Properties.RowNames = regionnames;

                if isfield(airway_res,'Wall_pct__3_8_')
                    T = addTableVarVal(T,'WallPct_3_8','WholeLung',airway_res.Wall_pct__3_8_);
                end

                fld = {'Wall_pct', '%0.2f' , {'Wall_pct',                     'WholeLung';...
                                              'Wall_pct_RightUpperLobe',      'RUL';...
                                              'Wall_pct_RightMidLobe',        'RML';...
                                              'Wall_pct_RightUpperLobePlus',  'RULplus';...
                                              'Wall_pct_RightLowerLobe',      'RLL';...
                                              'Wall_pct_LeftUpperLobe',       'LUL';...
                                              'Wall_pct_LeftUpperLobePlus',   'LULplus';...
                                              'Wall_pct_LeftLingula',         'LLi';...
                                              'Wall_pct_LeftLowerLobe',       'LLL'};...
                       'Pi10', '%0.2f' ,     {'Pi10',         'WholeLung';...
                                              'Pi10_RUL',     'RUL';...
                                              'Pi10_RML',     'RML';...
                                              'Pi10_RULplus', 'RULplus';...
                                              'Pi10_RLL',     'RLL';...
                                              'Pi10_LUL'      'LUL';...
                                              'Pi10_LULplus'  'RULplus';...
                                              'Pi10_LLi',     'LLi';...
                                              'Pi10_LLL',     'LLL'};...
                       'Pi15', '%0.2f' ,     {'Pi15',         'WholeLung'};...
                       'BEI',  '%0.2f' ,     {'BEI_Lung',     'WholeLung';...
                                              'BEI_Right',    'RL';...
                                              'BEI_Left',     'LL';...
                                              'BEI_RUL',      'RUL';...
                                              'BEI_RML',      'RML';...
                                              'BEI_RULplus',  'RULplus';...
                                              'BEI_RLL',      'RLL';...
                                              'BEI_LUL',      'LUL';...
                                              'BEI_LULplus',  'LULplus';...
                                              'BEI_LLi',      'LLi';...
                                              'BEI_LLL',      'LLL'};...
                       'BEI_gen', '%0.0f' ,  {'x1_Level_Lung',    'WholeLung';...
                                              'x1_Level_Right',   'RL';...
                                              'x1_Level_Left',    'LL';...
                                              'x1_Level_RUL',     'RUL';...
                                              'x1_Level_RML',     'RML';...
                                              'x1_Level_RULplus', 'RULplus';...
                                              'x1_Level_RLL',     'RLL';...
                                              'x1_Level_LUL',     'LUL';...
                                              'x1_Level_LLi',     'LLi';...
                                              'x1_Level_LULplus', 'LULplus';...
                                              'x1_Level_LLL',     'LLL'}};
                for icol = 1:size(fld,1)
                    for ifld = 1:size(fld{icol,3},1)
                        if isfield(airway_res,fld{icol,3}{ifld,1})
                            T = addTableVarVal(T,fld{icol,1},fld{icol,3}{ifld,2},airway_res.(fld{icol,3}{ifld,1}));
                        end
                    end
                end

                genstr = {'WT'};
                segstr = { 'seg',   4   ;...
                           'subseg',5:7 };
                lobestr =  {'Lung','Right','Left','RUL','RML','RUL+','RLL','LUL','LUL+','LLi','LLL'};
                lobestr2 = {'WholeLung','RL','LL','RUL','RML','RULplus','RLL','LUL','LULplus','LLi','LLL'};
                for i = 1:numel(genstr)
                    if isfield(airway_res,genstr{i})
                        for j = 1:size(segstr,1)
                            for k = 1:numel(lobestr)
                                if any(ismember(airway_res.(genstr{i}).Properties.VariableNames,lobestr{k}))
                                    segind = segstr{j,2};
                                    vals = airway_res.(genstr{i}).(lobestr{k});
                                    vals = vals(segind(segind<=numel(vals)));
                                    T = addTableVarVal(T,[genstr{i},'_',segstr{j,1}],lobestr2{k},mean(vals(vals>0)));
                                end
                            end
                        end
                    end
                end

            end
        end

        % Add table to results
        res = addTableVarVal(res,T);
    end

    % ScatterNet for AT on Exp CT scan
    try
        if opts.scatnetAT && img(1).flag
            writeLog(fn_log,'Air trapping map ... ');
            if exist(opts.fn.scatnetAT,'file')
                writeLog(fn_log,'from file\n');
                atMap = logical(cmi_load(1,img(1).info.d(1:3),opts.fn.scatnetAT));
            else
                writeLog(fn_log,'generating with ScatNet\n');
                atMap = ScatterNetAT(img(1).mat,logical(img(1).label),0);
                cmi_save(0,atMap,'ScatNet',img(1).info.fov,img(1).info.orient,opts.fn.scatnetAT);
            end
            T = CTlung_LobeStats(img(1).label,'scatnetAT','pct',atMap);
            % T = lobeLoop(img(1).label,@(mask,SN,img,str)tabulateScatNet(mask,SN,img,str),atMap,img(1).mat,'scatnetAT');
            res = addTableVarVal(res,T);
            clear atMap
        end
    catch err
        writeLog(fn_log,'ScatNet FAILED:\n%s\n',getReport(err));
    end
    
    % ScatNet for Emph on Ins CT scan
    try
        if opts.scatnetEmph && img(2).flag
            writeLog(fn_log,'scatnetEmph map ... ');
            if exist(opts.fn.scatnetEmph,'file')
                writeLog(fn_log,'from file\n');
                SNemph = logical(cmi_load(1,img(2).info.d(1:3),opts.fn.scatnetEmph));
            else
                writeLog(fn_log,'generating with ScatNet\n');
                SNemph = ScatterNetEMPH(img(2).mat,logical(img(2).label),0);
                cmi_save(0,SNemph,'ScatNet',img(2).info.fov,img(2).info.orient,opts.fn.scatnetEmph);
            end
            T = CTlung_LobeStats(img(2).label,'scatnetEmph','pct',SNemph);
            % T = lobeLoop(img(2).label,@(mask,SN,img,str)tabulateScatNet(mask,SN,img,str),SNemph,img(2).mat,'scatnetEmph');
            res = addTableVarVal(res,T);
            clear SNemph
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
        T = CTlung_Unreg('exp',img(1).mat,img(1).info.voxvol,img(1).label);
        res = addTableVarVal(res,T);
    end
    if opts.unreg && img(2).flag
        T = CTlung_Unreg('ins',img(2).mat,img(2).info.voxvol,img(2).label);
        res = addTableVarVal(res,T);
    end

    % Register I2E
    if img(1).flag && img(2).flag && (opts.reg || opts.prm || opts.tprm || opts.jac)
        ins_reg = [];
        elxdir = fullfile(procdir,sprintf('elastix_%s',ID));
        if exist(opts.fn.reg,'file')
            writeLog(fn_log,'Loading registered INS from file...\n');
            ins_reg = readNIFTI(opts.fn.reg);
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
        opts.fn.regMontage = fullfile(procdir,sprintf('%s_Reg_Montage.tif',res.ID{1}));
        QCmontage('reg',ins_reg(:,:,img(1).QCind),...
                        logical(img(1).label(:,:,img(1).QCind)),...
                        img(1).info.voxsz,...
                        opts.fn.regMontage,...
                        fullfile(opts.save_path,'Montage_reg.gif'));
        img(2) = [];

        % Jacobian analysis
        opts.fn.jac = fullfile(procdir,sprintf('%s.jac%s',res.ID{1},fn_ext));
        if opts.jac
            jac = [];
            if exist(opts.fn.jac,'file')
                writeLog(fn_log,'Analyzing Jacobian map...\n');
                jac = cmi_load(1,[],opts.fn.jac);
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
                        cmi_save(0,jac,'spatialJacobian',img(1).info.fov,img(1).info.orient,opts.fn.jac);
                        delete(fn_jac_elx);
                    else
                        writeLog(fn_log,'Jacobian processing failed.\n');
                    end
                end
            end
            if ~isempty(jac)
                T = CTlung_LobeStats(img.label,'Jac','mean',[],jac);
                % T = lobeLoop(img.label,@(mask,A,str)tabulateStats(mask,A,str),jac,'Jac');
                res = addTableVarVal(res,T);
            end
        end
        
        % ScatNet for Emph
        try
            if opts.scatnetEmph && img(1).flag
                writeLog(fn_log,'scatnetEmph map ... ');
                if exist(opts.fn.scatnetEmphReg,'file')
                    writeLog(fn_log,'from file\n');
                    SNemph = logical(cmi_load(1,img(1).info.d(1:3),opts.fn.scatnetEmphReg));
                else
                    writeLog(fn_log,'generating with ScatNet\n');
                    SNemph = ScatterNetEMPH(ins_reg,logical(img(1).label),0);
                    cmi_save(0,SNemph,'ScatNet',img(1).info.fov,img(1).info.orient,opts.fn.scatnetEmphReg);
                end
                T = CTlung_LobeStats(img.label,'scatnetEmphReg','pct',SNemph);
                % T = lobeLoop(img(1).label,@(mask,SN,img,str)tabulateScatNet(mask,SN,img,str),SNemph,ins_reg,'scatnetEmph');
                res = addTableVarVal(res,T);
                clear SNemph
            end
        catch err
            writeLog(fn_log,'ScatNet FAILED:\n%s\n',getReport(err));
        end
        
        % Blood density change map:
        if opts.dBlood
            dBlood = [];
            if exist(opts.fn.dBlood,'file')
                writeLog(fn_log,'Loading dBlood from file...\n');
                dBlood = readNIFTI(opts.fn.dBlood);
            elseif exist('ins_reg','var') && exist('jac','var')
                writeLog(fn_log,'Calculating dBlood...\n');
                dBlood = pipeline_blood_density(img(1).mat,ins_reg,jac,img(1).label);
                cmi_save(0,dBlood,{'dBlood'},img(1).info.fov,img(1).info.orient,opts.fn.dBlood);
            end
            if ~isempty(dBlood)
                T = CTlung_LobeStats(img.label,'dBlood','mean',[],dBlood);
                % T = lobeLoop(img.label,@(mask,A,str)tabulateStats(mask,A,str),dBlood,'dBlood');
                res = addTableVarVal(res,T);
            end
        end
        clear jac dBlood
        
        % PRM calculation
        prm10 = [];
        if opts.prm
            if exist(opts.fn.prm,'file')
                writeLog(fn_log,'Loading PRM from file...\n');
                prm10 = int8(readNIFTI(opts.fn.prm));
            else
                writeLog(fn_log,'Calculating PRM...\n');
                [prm10,~] = pipeline_PRM(img(1).mat,img(1).info,logical(img(1).label),ins_reg,...
                    opts.fn.prmScatter,fullfile(opts.save_path,'Montage_PRM_Scatter.gif'));

                % Save PRM
                writeLog(fn_log,'Saving PRM as NIFTI ... ');
                stat = cmi_save(0,prm10,{'PRM'},img.info.fov,img.info.orient,opts.fn.prm);
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
            T = CTlung_LobeStats(img.label,'PRM','pct',categorical(prm10,1:10,string(1:10)));
            % T = lobeLoop(img.label,@(mask,prm,flag)tabulatePRM(mask,prm,flag),prm10,1);
            res = addTableVarVal(res,T);

            % map full PRM values (1:10) to (norm,fsad,emph,pd,ns)
            prmtranslate = {'Norm', 'fSAD', 'Emph', 'PD',       'NS';...
                            [1,2],  3,      [4,5],  [8,9,10],   6    };
            prm5 = int8(zeros(size(prm10)));
            for i = 1:size(prmtranslate,2)
                prm5(ismember(prm10,prmtranslate{2,i})) = i;
            end
            clear prm10

            % QC PRM
            writeLog(fn_log,'Generating PRM Montage ...\n');
            QCmontage('prm',img.mat(:,:,img(1).QCind),...
                            double(prm5(:,:,img(1).QCind)),...
                            img.info.voxsz,...
                            opts.fn.prmMontage,...
                            fullfile(opts.save_path,'Montage_PRM.gif'));

            % Tabulate 5-color PRM results
            writeLog(fn_log,'Tabulating 5-color PRM results...\n');
            % T = lobeLoop(img.label,@(mask,prm,flag)tabulatePRM(mask,prm,flag),prm5,0);
            T = CTlung_LobeStats(img.label,'PRM','pct',categorical(prm5,1:5,string(prmlabel)));
            res = addTableVarVal(res,T);

            % Calculate tPRM
            if opts.tprm
                % Check if tPRM maps already exist
                fn_chk = all(cellfun(@(x)exist(x,'file'),opts.fn.tprm),"all");

                t = tic;
                if ~fn_chk
                    % Need to calculate gridded MF
                    writeLog(fn_log,'Generating tPRM ...\n');
                    p = minkowskiFun(prm5,'thresh',1:4,...
                        'tmode','==',...
                        'n',10*ones(1,3),...
                        'gridsp',5,...
                        'voxsz',img.info.voxsz,...
                        'mask',logical(img.label),...
                        'prog',0);
                else
                    writeLog(fn_log,'Loading tPRM from files ...\n');
                end
                clear prm5;

                % Analyze tPRM means
                for iprm = 1:4
                    for imf = 1:numel(mflabel)
                        str = prmlabel(iprm) + '_' + upper(mflabel(imf));
                        writeLog(fn_log,'   %s ... ',str);
                        if fn_chk
                            writeLog(fn_log,'from file ... ');
                            tprm = readNIFTI(opts.fn.tprm(imf,iprm));
                        else
                            writeLog(fn_log,'interpolating ... ');
                            tprm = grid2img(p.MF(iprm,imf,:),p.ind,p.mask,3,1);
                            writeLog(fn_log,'Saving NIFTI ... ');
                            cmi_save(0,single(tprm),{char(str)},img.info.fov,img.info.orient,opts.fn.tprm{imf,iprm});
                        end
                        writeLog(fn_log,'Tabulating means\n');
                        T = CTlung_LobeStats(img.label, "tPRM_" + str, 'mean', [], tprm);
                        res = addTableVarVal(res,T);
                    end
                end
                writeLog(fn_log,'... tPRM complete (%s)\n',datetime([0,0,0,0,0,toc(t)],'Format','HH:mm:ss'));
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
    
writeLog(fn_log,'Pipeline total time = %s\n',datetime([0,0,0,0,0,toc(tt)],'Format','HH:mm:ss'))

function T = tabulateStats(mask,A,str)
    vname = strcat(str,{'_mean','_var'});
    nv = numel(vname);
    T = table('Size',[1,nv],'VariableTypes',repmat({'double'},1,nv),'VariableNames',vname);
    T.(vname{1}) = mean(A(mask));
    T.(vname{2}) = var(A(mask));
    
function T = tabulateScatNet(mask,scatnet,img,str)
    vname = strcat(str,{'_pct','_mean'});
    nv = numel(vname);
    T = table('Size',[1,nv],'VariableTypes',repmat({'double'},1,nv),'VariableNames',vname);
    T.(vname{1}) = nnz(scatnet(mask))/nnz(mask)*100;
    T.(vname{2}) = mean(img(mask & scatnet));

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
    