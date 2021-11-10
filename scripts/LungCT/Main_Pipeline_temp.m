function T = Main_Pipeline_temp(catalog)

% Input: catalog = cell array of dicom information in a directory generated
%                   by the "dcmCatalog()" script

% Pipeline consists of:
% 1. Select data for processing based on catalogDICOM script .csv file
% 2. Identify image key information and read images
% 3. Generate lung segmentations
% 4. Identify Exp/Ins images and save as .nii.gz
% 5. Quantify unregistered lung statistics
% 2. Enqueue elastix registration

%% Determine data for processing:
% cases is structure array:
%   cases(i).PatientID
%           .StudyDate
%           .Scans - Structure containing scan info for selected data
    if nargin==0
        catalog = [];
    end
    [cases, home_pwd] = catalog_select(catalog);

    %% Initialize results table:
    N = length(cases);
    vars = {'key',          'int',...
            'ID',           'cellstr';...
            'Exp_DICOM',    'cellstr';...
            'Ins_DICOM',    'cellstr'...
            'Proc_Location','cellstr'};
    T = table('Size',[N,size(vars,1)],'VariableTypes',vars(:,2),'VariableNames',vars(:,1));

    %% Loop over cases
    fn_ext = '.nii.gz';
    h1 = waitbar(0, 'Analyze Exp and Ins data');
    for i = 1:N
        ID = sprintf('%s_%s',cases(i).PatientName,cases(i).StudyDate);
        waitbar(i/N,h1,[num2str(round(100*(i-1)/N,1)),'% Complete: Load Data and Masks for ',ID])

        %% Establish relevant filenmes:
        procdir = fullfile(home_pwd,'ProcessedData',ID);
        dcmexp = cases(i).Scans(strcmp({cases(i).Scans(:).Tag},'Exp')).Directory;
        dcmins = cases(i).Scans(strcmp({cases(i).Scans(:).Tag},'Ins')).Directory;
        fn.exp.name = sprintf('%s.exp',ID);
        fn.ins.name = sprintf('%s.ins',ID);
        fn.explabel.name = sprintf('%s.exp.label',ID);
        fn.inslabel.name = sprintf('%s.ins.label',ID);
        fn.scatnet.name = sprintf('%s.scatnet',ID);
        fn.reg.name = sprintf('%s.ins.reg',ID);
        fn.prm.name = sprintf('%s.prm',ID);
        fn.tprm.name = cellfun(@(x)sprintf('%s.tprm.%s',ID,x),...
            append({'norm','fsad','emph','pd'},'.',{'v';'s';'b';'x'}),'UniformOutput',false);
        fld = fieldnames(fn);
        for ifn = 1:numel(fld)
            fn.(fld{ifn}).path = fullfile(procdir,strcat(fn.(fld{ifn}).name,fn_ext));
        end
        elxdir = fullfile(procdir,sprintf('elxreg_%s',fname_Ins));
        if ~isfolder(elxdir)
            mkdir(elxdir);
        end

        T.key = i;
        T.ID{i} = ID;
        T.Proc_Location{i} = procdir;
        T.Exp_DICOM{i} = dcmexp;
        T.Ins_DICOM{i} = dcmins;

        %% Initialize image struct
        img = struct('mat',{[],[]},'info',{[],[]},'label',{[],[]});

        %% Read selected DICOM data:
        fprintf('\nReading image data from file ... ID = %s\n',ID)
        if exist(fn.exp.path,'file') && exist(fn.ins.path,'file')
            fprintf('   ... from NiFTi\n');

            fprintf('   ... Reading Exp\n');
            [img(1).mat,label,fov,orient,info] = readNIFTI(fn.exp.path);
            d = size(img(1).mat);
            img(1).info = struct('label',label,'fov',fov,'orient',orient,'d',d,'voxsz',fov./d,'natinfo',info);
            img(1).info.voxvol = prod(img(1).info.voxsz);

            fprintf('   ... Reading Ins\n');
            [img(2).mat,label,fov,orient,info] = readNIFTI(fn.ins.path);
            d = size(img(2).mat);
            img(2).info = struct('label',label,'fov',fov,'orient',orient,'d',d,'voxsz',fov./d,'natinfo',info);
            img(2).info.voxvol = prod(img(2).info.voxsz);        

            check_EI = false;
        else
            fprintf('   ... from DICOM\n');

            [img(1).mat,label,fov,orient,info] = readDICOM(fn.exp.path);
            d = size(img(1).mat);
            img(1).info = struct('label',label,'fov',fov,'orient',orient,'d',d,'voxsz',fov./d,'natinfo',info);

            [img(2).mat,label,fov,orient,info] = readDICOM(fn.ins.path);
            d = size(img(2).mat);
            img(2).info = struct('label',label,'fov',fov,'orient',orient,'d',d,'voxsz',fov./d,'natinfo',info);

            check_EI = true;
        end

        %% Generate Lung Segmentation [This is when VOI don't exist]
        if exist(fn.explabel.path,'file') && exist(fn.inslabel.path,'file')
            fprintf('   Reading VOIs from file\n')
            img(1).label = readNIFTI(fn.explabel.path);
            img(2).label = readNIFTI(fn.inslabel.path);
        else
            fprintf('   Generating VOIs\n')
            img(1).label = DL_lung_segmetation(img(1).mat);
            img(2).label = DL_lung_segmetation(img(2).mat);
        end

        %% Identify Exp and Ins using lung volume; used for determining file name
        img(1).info.voxvol = prod(img(1).info.voxsz);
        img(2).info.voxvol = prod(img(2).info.voxsz);
        if check_EI
            %   ** Need to double-check in case of mislabel
            if (nnz(img(1).label)*img(1).info.voxvol) > (nnz(img(2).label)*img(2).info.voxvol)
                fprintf('Swapping INS/EXP due to lung volume\n');
                img = img([2,1]);
                T.Exp_Location{i} = dcmins;
                T.Ins_Location{i} = dcmexp;
            end
            %% Save nii.gz files using ID and Tag
            saveNIFTI(fn.exp.path,img(1).mat,img(1).info.label,img(1).info.fov,img(1).info.orient);
            saveNIFTI(fn.ins.path,img(2).mat,img(2).info.label,img(2).info.fov,img(2).info.orient);
            saveNIFTI(fn.explabel.path,img(1).label,img(1).info.label,img(1).info.fov,img(1).info.orient);
            saveNIFTI(fn.inslabel.path,img(2).label,img(2).info.label,img(2).info.fov,img(2).info.orient);
        end

        %% QC segmentation
        ind = 10:10:img(1).info.d(3);
        QCmontage('seg',cat(4,img(1).mat(:,:,ind),img(1).label(:,:,ind)),img(1).info.voxsz,...
            fullfile(procdir,sprintf('%s_Montage',fname_Exp)));
        ind = 10:10:img(2).info.d(3);
        QCmontage('seg',cat(4,img(2).mat(:,:,ind),img(2).label(:,:,ind)),img(2).info.voxsz,...
            fullfile(procdir,sprintf('%s_Montage',fname_Ins)));

        %% ScatterNet for AT on Exp CT scan
        if exist(fn.scatnet.path,'file')
            atMap = cmi_load(1,img(1).info.d(1:3),fn.scatnet.path);
        else
            atMap = ScatNet(img(1).mat,logical(img(1).label),0);
            cmi_save(0,atMap,'ScatNet',img(1).info.fov,img(1).info.orient,fn.scatnet.path);
        end

        %% Quantify unregistered CT scans
        fprintf('\n   Quantifying unregistered statistics\n');
        S = CTlung_Unreg('exp',img(1).mat,img(1).info.voxvol,img(1).label,atMap);
        clear atMap;
        for ilab = 1:numel(S)
            tstr = sprintf('_%u',S(ilab).tag);
            T.(['Exp_Vol',tstr])(i) = S(ilab).vol;
            T.(['Exp_HU',tstr])(i) = S(ilab).mean;
            T.(['Exp_856',tstr])(i) = S(ilab).exp856;
            T.(['Exp_SNpct',tstr])(i) = S(ilab).SNpct;
            T.(['Exp_SNmean',tstr])(i) = S(ilab).SNmean;
        end
        S = CTlung_Unreg('ins',img(2).mat,img(2).info.voxvol,img(2).label);
        for ilab = 1:numel(S)
            tstr = sprintf('_%u',S(ilab).tag);
            T.(['Ins_Vol',tstr])(i) = S(ilab).vol;
            T.(['Ins_HU',tstr])(i) = S(ilab).mean;
            T.(['Ins_950',tstr])(i) = S(ilab).ins950;
            T.(['Ins_810',tstr])(i) = S(ilab).ins810;
        end

        %% Register I2E
        if exist(fn.reg.path,'file')
            fprintf('Loading registered INS from file...\n')
        else
            fprintf('Performing registration...\n')
            lungreg_BH( img(1).mat, img(1).info, logical(img(2).label),...
                        img(2).mat, img(2).info, logical(img(1).label),...
                        elxdir, ID, false);
            % Move registered file out to procdir
            movefile(fullfile(elxdir,sprintf('%s_R.nii',img(2).info.name)),fn.reg.path);
        end
        img(2) = [];
        fprintf('Loading ins.reg file...\n');
        ins_reg = readNIFTI(fn.reg.path);

        %% QC registration
        ind = 10:10:img(1).info.d(3);
        QCmontage('reg',cat(4,ins_reg(:,:,ind),img(1).label(:,:,ind)),img(1).info.voxsz,...
            fullfile(procdir,sprintf('%s_Montage',fn.reg.name)));

        %% PRM calculation
        if exist(fullfile(procdir,[fname_PRM,fn_ext]),'file')
            fprintf('Loading PRM from file...\n');
            prm = int8(readNIFTI(fn.prm.path));
        else
            fprintf('Calculating PRM...\n');
            [prm,prminfo] = pipeline_PRM(img(1).mat,img(1).info,logical(img(1).label),ins_reg);

            % QC PRM
            print(prminfo.fig,fullfile(procdir,sprintf('%s_PRM_Scatter',ID)),'-dtiff');
            QCmontage('prm',cat(4,img.mat,prm),img.info.voxsz,fullfile(procdir,sprintf('%s_PRM_Montage',ID)));

            % Save PRM
            stat = cmi_save(0,prm,{fn.prm.name},img.info.fov,img.info.orient,fn.prm.path);
            if stat
                fprintf('  PRM saved\n');
            else
                fprintf('  Could not save PRM to file.\n');
            end
        end
        clear ins_reg;
        % map full PRM values (1:10) to (norm,fsad,emph,pd,ns)
        prmlabel = {'Norm','fSAD','Emph','PD','NS'};
        prm(ismember(prm,1:2)) = 1; % Norm
        prm(prm==3)            = 2; % fSAD
        prm(ismember(prm,4:5)) = 3; % Emph
        prm(prm==6)            = 5; % NS
        prm(ismember(prm,8:10))= 4; % PD

        %% Tabulate PRM results:
        BW = logical(img.label);
        ulab = unique(img.label(BW));
        nlab = numel(ulab);
        tstr = '';
        for ilab = 1:(nlab+(nlab>1))
            if ilab>1
                BW = img.label == (ulab(ilab-1));
                tstr = sprintf('_%u',ulab(ilab-1));
            end
            for iprm = 1:numel(prmlabel)
                T.(['PRM_',prmlabel{iprm,1},tstr])(i) = nnz(prm(BW)==iprm);
            end
        end

        %% Calculate tPRM
        if all(cellfun(@(x)exist(x,'file'),fn.tprm.path))
            fprintf('Loading tPRM from file...\n');
            clear prm;
        else
            t = tic;
            fprintf('Generating tPRM...\n');
            
            %% Calculate MF values
            p = minkowskiFun(prm,'thresh',1:4,...
                                 'tmode','==',...
                                 'n',10*ones(1,3),...
                                 'gridsp',5,...
                                 'voxsz',img.info.voxsz,...
                                 'mask',logical(img.label),...
                                 'prog',1);
            clear prm;
            
            %% Interpolate to maps
            for ithresh = 1:size(p.MF,1)
                for imf = 1:size(p.MF,2)
                    str_th = prmlabel{ithresh};
                    str_mf = p.labels{imf};
                    
                    % Interpolate to image space
                    tstr = [str_th,'_',str_mf];
                    fprintf('   Interpolating %s\n',tstr);
                    tprm = grid2img(p.MF(),p.ind,p.mask,3,1);
                    
                    % Save tPRM image
                    fname = join({ID,'tprm',lower(str_th),lower(str_mf)},'.');
                    fpath = fullfile(procdir,[fname,fn_ext]);
                    fprintf('   Saving: %s\n',fname);
                    cmi_save(0,tprm,{fname},img.info.fov,img.info.orient,fpath);
                    
                    % Tabulate statistics
                    BW = logical(img.label);
                    fstr = '';
                    for ilab = 1:(nlab+(nlab>1))
                        if ilab>1
                            BW = img.label==ulab(ilab-1);
                            fstr = sprintf('_%u',ulab(ilab-1));
                        end
                        T.([tstr,fstr])(i) = mean(tprm(BW));
                    end
                end
            end
            
            fprintf('... tPRM complete (%s)\n',datestr(duration(0,0,toc(t)),'HH:MM:SS'));
            
            fprintf('Saving tPRM results...\n');
            for itprm = 1:size(tprm,4)
            end
        end


    end
    delete(h1);
    
    %% Save Results Table:
    writetable(T,fullfile(procdir,fullfile(home_pwd,'ProcessedData','Results.csv')))
    
end

function QCmontage(tag,mat,voxsz,fname)
    nf = size(mat,3);

    % Initialize figure:
    hf = figure('visible','off');
    ha = axes(hf,'DataAspectRatio',voxsz);
    im = imagesc(mat(:,:,1,1),clim);
    axis image
    set(ha,'Units','pixels');
    npos = plotboxpos(ha);
    npos(1:2) = npos(1:2)+1;
    npos(3:4) = npos(3:4)-2;
    F = getframe(hf,npos);
    montsz = [round(sqrt(nf)),nan];

    % Generate frames:
    switch tag
        case {'seg','reg'}
            colormap(hf,'gray');
            hroi = plot(ha,1,1,'*m','MarkerSize',10);
            for i = 1:nf
                im.CData = mat(:,:,i,1);
                [voir,voic] = find(edge(mat(:,:,i,2),'Canny'));
                set(hroi,'XData',voic,'YData',voir);
                F(i) = getframe(hf,npos);
            end
        case 'prm'
            colormap(hf,[ gray(256) ;...
                          0 1 0 ;... % Norm
                          1 1 0 ;... % fSAD
                          1 0 0 ;... % Emph
                          1 0 1 ;... % PD
                          1 1 1 ]);  % NS
            tprm = mat(:,:,:,2);
            mat = (max(min(mat(:,:,:,1),0),-1000)+1000)*256/1000 - 0.5;
            mat(tprm>0) = tprm + 256;
            ha.CLim = [0 261];
            for i = 1:nf
                im.CData = mat(:,:,i);
                F(i) = getframe(hf,npos);
            end
    end

    % Generate montage
    mat = zeros(size(F(1).cdata,1),size(F(1).cdata,2),3,nf);
    for i = 1:nf
        [im,map] = frame2im(F(i));
        if ~isempty(map)
            im = ind2rgb(im,map);
        end
        mat(:,:,:,i) = double(im);
    end
    figure(hf),montage(mat/255,'Size',montsz);

    %% Print the figure:
    print(hf,fname,'-dtiff');
    delete(hf);
end

