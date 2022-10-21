function res = Main_Pipeline_GL_sub(procdir,opts)
% Pipeline processing on Great Lakes

tt = tic;

try

    fn_ext = '.nii.gz';
    opts.fn_log = fullfile(procdir,'pipeline_log.txt');

    writeLog(opts.fn_log,'\n\nStarting Pipeline Process : %s\n',datestr(datetime('now'),0));
    
    % Load results file
    [~,ID] = fileparts(procdir);
    regionnames = {'WholeLung','RL','LL','RUL','RML','RULplus','RLL','LUL','LULplus','LLi','LLL'};
    Nr = numel(regionnames);
    fn_res = fullfile(fullfile(procdir,[ID,'_PipelineResults.csv']));
    if exist(fn_res,'file')
        res = readtable(fn_res,'Delimiter',',');
    else
        varnames = {'ID','Exp_DICOM','Ins_DICOM'};
        Nv = numel(varnames);
        res = table('Size',[Nr,Nv],'VariableTypes',repmat({'cellstr'},1,Nv),'VariableNames',varnames);
        res.ID(:) = {ID};
    end
    res.Properties.RowNames = regionnames;
    res.ROI = regionnames';

    % Load images
    img = struct('flag',{false,false},'mat',{[],[]},'info',{[],[]},'label',{[],[]});
    fn_exp = fullfile(procdir,[ID,'.exp',fn_ext]);
    if exist(fn_exp,'file')
        writeLog(opts.fn_log,'Reading EXP image ...\n');
        [img(1).mat,label,fov,orient,info] = cmi_load(1,[],fn_exp);
        d = size(img(1).mat);
        voxsz = fov./d;
        img(1).info = struct('label',label,'fov',fov,'orient',orient,'d',d,'voxsz',voxsz,'natinfo',info,...
            'voxvol',prod(voxsz),'name',label);
        img(1).flag = true;
    end
    fn_ins = fullfile(procdir,[ID,'.ins',fn_ext]);
    if exist(fn_ins,'file')
        writeLog(opts.fn_log,'Reading INS image ...\n');
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
            img(1).label = getRespiratoryOrgans(medfilt2_3(img(1).mat));
        end
    end
    fn_ins_seg = fullfile(procdir,[ID,'.ins.label',fn_ext]);
    if img(2).flag
        if exist(fn_ins_seg,'file')
            img(2).label = cmi_load(1,[],fn_ins_seg);
        else
            img(2).label = getRespiratoryOrgans(medfilt2_3(img(2).mat));
        end
    end

    % QC segmentation
    if img(1).flag
        writeLog(opts.fn_log,'Generating EXP montage...\n');
        ind = 10:10:img(1).info.d(3);
        cdata = QCmontage('seg',cat(4,img(1).mat(:,:,ind),logical(img(1).label(:,:,ind))),img(1).info.voxsz,...
            fullfile(procdir,sprintf('%s_Montage',img(1).info.label)));
    end
    if img(2).flag
        writeLog(opts.fn_log,'Generating INSP montage...\n');
        ind = 10:10:img(2).info.d(3);
        cdata = QCmontage('seg',cat(4,img(2).mat(:,:,ind),logical(img(2).label(:,:,ind))),img(2).info.voxsz,...
            fullfile(procdir,sprintf('%s_Montage',img(2).info.label)));
    end

    % Airways
    if opts.airway
        ei_str = {'Exp','Ins'};
        for itag = 1:2
            if img(itag).flag
                ydir = fullfile(procdir,['yacta_',ID,'_',ei_str{itag}]);
                writeLog(opts.fn_log,'YACTA directory: %s\n',ydir);
                airway_res = readYACTAairways(ydir);
                if ~isempty(airway_res)

                    res = addTableVarVal(res,'WallPct_3_8','WholeLung',airway_res.Wall_pct__3_8_);

                    res = addTableVarVal(res,'Wall_pct',{'WholeLung','RUL','RML','RULplus','RLL','LUL','LULplus','LLi','LLL'},...
                        [ airway_res.Wall_pct;...
                        airway_res.Wall_pct_RightUpperLobe;...
                        airway_res.Wall_pct_RightMidLobe;...
                        airway_res.Wall_pct_RightUpperLobePlus;...
                        airway_res.Wall_pct_RightLowerLobe;...
                        airway_res.Wall_pct_LeftUpperLobe;...
                        airway_res.Wall_pct_LeftUpperLobePlus;...
                        airway_res.Wall_pct_LeftLingula;...
                        airway_res.Wall_pct_LeftLowerLobe ]);

                    res = addTableVarVal(res,'Pi10',{'WholeLung','RUL','RML','RULplus','RLL','LUL','LULplus','LLi','LLL'},...
                        [ airway_res.Pi10;...
                        airway_res.Pi10_RUL;...
                        airway_res.Pi10_RML;...
                        airway_res.Pi10_RULplus;...
                        airway_res.Pi10_RLL;...
                        airway_res.Pi10_LUL;...
                        airway_res.Pi10_LULplus;...
                        airway_res.Pi10_LLi;...
                        airway_res.Pi10_LLL ]);

                    res = addTableVarVal(res,'Pi15','WholeLung',airway_res.Pi15);

                    genstr = {'WT'};
                    segstr = { 'seg',   4   ;...
                        'subseg',5:7 };
                    lobestr =  {'Right','Left','RUL','RML','RUL+','RLL','LUL','LUL+','LLi','LLL'};
                    lobestr2 = {'RL','LL','RUL','RML','RULplus','RLL','LUL','LULplus','LLi','LLL'};
                    for i = 1:numel(genstr)
                        for j = 1:2
                            for k = 1:numel(lobestr)
                                segind = segstr{j,2};
                                vals = airway_res.(genstr{i}).(lobestr{k});
                                vals = vals(segind(segind<=numel(vals)));
                                res = addTableVarVal(res,[genstr{i},'_',segstr{j,1}],lobestr2{k},mean(vals(vals>0)));
                            end
                        end
                    end

                    res = addTableVarVal(res,'BEI',{'WholeLung','RL','LL','RUL','RML','RULplus','RLL','LUL','LULplus','LLi','LLL'},...
                        [ airway_res.BEI_Lung;...
                        airway_res.BEI_Right;...
                        airway_res.BEI_Left;...
                        airway_res.BEI_RUL;...
                        airway_res.BEI_RML;...
                        airway_res.BEI_RULplus;...
                        airway_res.BEI_RLL;...
                        airway_res.BEI_LUL;...
                        airway_res.BEI_LULplus;...
                        airway_res.BEI_LLi;...
                        airway_res.BEI_LLL ]);
                end
            end
        end
    end

    % ScatterNet for AT on Exp CT scan
    if opts.scatnet && img(1).flag
        fn_scatnet = fullfile(procdir,sprintf('%s.%s%s',res.ID{1},'scatnet',fn_ext));
        writeLog(opts.fn_log,'Air trapping map ... ');
        if exist(fn_scatnet,'file')
            writeLog(opts.fn_log,'from file\n');
            atMap = cmi_load(1,img(1).info.d(1:3),fn_scatnet);
        else
            writeLog(opts.fn_log,'generating with ScatNet\n');
            atMap = ScatNet(img(1).mat,logical(img(1).label),0);
            cmi_save(0,atMap,'ScatNet',img(1).info.fov,img(1).info.orient,fn_scatnet);
        end
    end

    % Vessel analysis
    if opts.vessel && img(2).flag
        writeLog(opts.fn_log,'Vessel analysis ...\n');
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
        T = pipeline_vesselseg( tins, tseg, tinfo, procdir, opts, opts.fn_log);
        res = addTableVarVal(res,T);
    end

    % Quantify unregistered CT scans
    writeLog(opts.fn_log,'Quantifying unregistered statistics\n');
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
        if exist(fn_reg,'file')
            writeLog(opts.fn_log,'Loading registered INS from file...\n');
            ins_reg = readNIFTI(fn_reg);
        elseif opts.reg
            elxdir = fullfile(procdir,sprintf('elastix_%s',ID));
            writeLog(opts.fn_log,'Performing registration...\nelxdir : %s\n',elxdir)
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
                    writeLog(opts.fn_log,'Re-saving: %s\n',fn_move{i,1});
                    timg = cmi_load(1,[],fullfile(elxdir,fn(1).name));
                    cmi_save(0,timg,fn_move{i,1},img(1).info.fov,img(1).info.orient,...
                        fullfile(procdir,sprintf('%s%s%s',res.ID{1},fn_move{i,2},fn_ext)));
                    if i==nf
                        ins_reg = timg;
                        clear timg;
                    end
                end
            end
            
        end

        % QC registration
        writeLog(opts.fn_log,'Saving Registration Montage ...\n');
        ind = 10:10:img(1).info.d(3);
        QCmontage('reg',cat(4,ins_reg(:,:,ind),img(1).label(:,:,ind)),img(1).info.voxsz,...
            fullfile(procdir,sprintf('%s_Reg_Montage',res.ID{1})));
        img(2) = [];

        fn_jac = fullfile(procdir,sprintf('%s.jac%s',res.ID{1},fn_ext));
        if opts.jac
            if exist(fn_jac,'file')
                writeLog(opts.fn_log,'Analyzing Jacobian map...\n');
                jac = cmi_load(1,[],fn_jac);
                T = lobeLoop(img.label,@(mask,jac)tabulateJac(mask,jac),jac);
                res = addTableVarVal(res,T);
            else
                writeLog(opts.fn_log,'Jacobian file not found.\n');
            end
        end
        
        % PRM calculation
        prm10 = [];
        fn_PRM = fullfile(procdir,sprintf('%s.%s%s',res.ID{1},'prm',fn_ext));
        if exist(fn_PRM,'file')
            writeLog(opts.fn_log,'Loading PRM from file...\n');
            prm10 = int8(readNIFTI(fn_PRM));
        elseif opts.prm
            writeLog(opts.fn_log,'Calculating PRM...\n');
            [prm10,~] = pipeline_PRM(img(1).mat,img(1).info,logical(img(1).label),ins_reg,...
                fullfile(procdir,sprintf('%s_PRM_Scatter',res.ID{1})));

            % Save PRM
            writeLog(opts.fn_log,'Saving PRM as NIFTI ... ');
            stat = cmi_save(0,prm10,{'PRM'},img.info.fov,img.info.orient,fn_PRM);
            if stat
                writeLog(opts.fn_log,'  PRM saved\n');
            else
                writeLog(opts.fn_log,'  Could not save PRM to file.\n');
            end
        end
        clear ins_reg;

        if ~isempty(prm10)
            % Tabulate 10-color PRM results
            writeLog(opts.fn_log,'Tabulating 10-color PRM results...\n');
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
            writeLog(opts.fn_log,'Generating PRM Montage ...\n');
            QCmontage('prm',cat(4,img.mat(:,:,ind),double(prm5(:,:,ind))),...
                img.info.voxsz,fullfile(procdir,sprintf('%s_PRM_Montage',res.ID{1})));

            % Tabulate 5-color PRM results
            writeLog(opts.fn_log,'Tabulating 5-color PRM results...\n');
            T = lobeLoop(img.label,@(mask,prm,flag)tabulatePRM(mask,prm,flag),prm5,0);
            res = addTableVarVal(res,T);

            % Calculate tPRM
            prmlabel = ["norm","fsad","emph","pd"];
            mflabel = ["v","s","b","x"];
            fn_tprm = fullfile(procdir,...
                string(res.ID{1})+".tprm."+prmlabel+"."+mflabel'+string(fn_ext));
            if all(cellfun(@(x)exist(x,'file'),fn_tprm))
                writeLog(opts.fn_log,'Loading tPRM from files ...\n');
                clear prm5;
                for iprm = 1:numel(prmlabel)
                    for imf = 1:numel(mflabel)
                        writeLog(opts.fn_log,'   %s - %s\n',prmlabel(iprm),mflabel(imf))
                        tprm = readNIFTI(fn_tprm(imf,iprm));
                        str = prmlabel(iprm)+'_'+upper(mflabel(imf));
                        T = lobeLoop(img.label,@(mask,tprm,str)tabulateTPRM(mask,tprm,str),tprm,str);
                        res = addTableVarVal(res,T);
                    end
                end
            elseif opts.tprm
                t = tic;
                writeLog(opts.fn_log,'Generating tPRM ...\n');

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
                        writeLog(opts.fn_log,'   %s - %s\n',prmlabel(ithresh),mflabel(imf));

                        % Interpolate to image space
                        tstr = prmlabel(ithresh) + '.' + mflabel(imf);
                        writeLog(opts.fn_log,'       Interpolating\n');
                        tprm = grid2img(p.MF(ithresh,imf,:),p.ind,p.mask,3,1);

                        % Save tPRM image
                        writeLog(opts.fn_log,'       Saving NIFTI\n');
                        cmi_save(0,single(tprm),{char(tstr)},img.info.fov,img.info.orient,char(fn_tprm(imf,ithresh)));

                        % Tabulate statistics
                        writeLog(opts.fn_log,'       Tabulating means\n');
                        str = prmlabel(ithresh)+'_'+upper(mflabel(imf));
                        T = lobeLoop(img.label,@(mask,tprm,str)tabulateTPRM(mask,tprm,str),tprm,str);
                        res = addTableVarVal(res,T);
                    end
                end
                writeLog(opts.fn_log,'... tPRM complete (%s)\n',datestr(duration(0,0,toc(t)),'HH:MM:SS'));

            end
        end
    end

    % Save Results Table:
    if istable(res)
        writetable(res,fullfile(procdir,[res.ID{1},'_PipelineResults.csv']));
    end
    
catch err
    writeLog(opts.fn_log,'Pipeline ERROR:\n%s',getReport(err));
end

writeLog(opts.fn_log,'Pipeline total time = %s\n',datestr(duration(0,0,toc(tt)),'HH:MM:SS'))

function T = tabulateJac(mask,jac)
    vname = {'Jac_mean'};
    nv = numel(vname);
    T = table('Size',[1,nv],'VariableTypes',repmat({'double'},1,nv),'VariableNames',vname);
    T.Jac_mean = mean(jac(mask));

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
    
function T = addTableVarVal(T,varargin)
    if nargin==2 && istable(varargin{1}) && strcmp(varargin{1}.Properties.VariableNames{1},'LOBE')
        vals = varargin{1};
        regionstr = vals.LOBE;
        vals = removevars(vals,'LOBE');
        varstr = vals.Properties.VariableNames;
    elseif nargin==4
        varstr = varargin{1};
        regionstr = varargin{2};
        vals = varargin{3};
        if ischar(regionstr)
            regionstr = {regionstr};
        end
        if ischar(varstr)
            varstr = {varstr};
        end
        if numel(varstr)<size(vals,2)
            varstr = strcat(varstr,'_',cellfun(@num2str,num2cell(1:size(vals,2)),'UniformOutput',false));
        end
    else
        return;
    end
    if numel(regionstr)==size(vals,1)
        for i = 1:numel(varstr)
            if istable(vals)
                tvals = vals.(varstr{i});
            else
                tvals = vals(:,i);
            end
            if iscell(tvals)
                defval = {''};
            elseif isnumeric(tvals)
                defval = nan;
            else
                defval = {[]};
            end
            if ~ismember(varstr{i},T.Properties.VariableNames)
                T = addvars(T,repmat(defval,size(T,1),1),'NewVariableNames',varstr(i));
            end
            T.(varstr{i})(regionstr) = tvals;
        end
    end
    