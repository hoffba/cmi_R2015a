function CTlung_Pipeline_local(ID,expfname,insfname,procdir,opts)
% Local execution for YACTA segmentation and airways analysis
% * procdir must be on Turbo for Great Lakes to have access to it

t = tic;
try
    
    % Save data sources to table in processing directory for future reference
    T = table('Size',[2,2],'VariableTypes',{'cellstr','cellstr'},'VariableNames',{'Tag','DataPath'});

    % Initialize folders and log file
    if ~isfolder(procdir)
        mkdir(procdir);
    end
    fn_log = fullfile(procdir,sprintf('pipeline_log.txt'));
    writeLog(fn_log,'%s\n',ID);

    % Initialize parameters and results struct
    fn_ext = '.nii.gz';
    img = struct('flag',[false,false],'mat',{[],[]},'info',{[],[]},'label',{'',''},'sourcepath',{'Unknown','Unknown'});

    fn = fullfile(procdir,string(ID)+[".exp",".exp.label";".ins",".ins.label"]+fn_ext);
    fnflag = cellfun(@(x)exist(x,'file'),fn);
    if any(~fnflag(:,1))
        % If either image is missing, must redo all from DICOM
        fnflag(:) = false;
    end

    % Load CT data
    tagstr = {'Exp','Ins'};
    fn_in = {expfname,insfname};
    for i = 1:2
        orientchk = false;
        % CT from origin files
        writeLog(fn_log,'%s : CT ... ',tagstr{i});
        if fnflag(i,1)
            writeLog(fn_log,'file found: %s\n',fn{i,1});
            if ~fnflag(i,2)
                % Need to load image for segmentation
                writeLog(fn_log,'   Loading for segmentation.\n');
                [img(i).mat,label,fov,orient,info] = cmi_load(1,[],fn_in{i});
            end
        else
            orientchk = true;
            writeLog(fn_log,'   loading data ... ',fn{i,1});
            img(i).sourcepath = fn_in{i};
            if isfolder(fn_in{i})
                writeLog(fn_log,'from DICOM\n');
                [img(i).mat,label,fov,orient,info] = readDICOM(fn_in{i},[],'noprompt');
            elseif ~isempty(fn_in{i})
                writeLog(fn_log,'from file\n')
                [img(i).mat,label,fov,orient,info] = cmi_load(1,[],fn_in{i});
            else
                writeLog(fn_log,'image not loaded.\n');
            end
        end

        img(i).flag = ~isempty(img(i).mat);
        if img(i).flag
            % Set image info
            d = size(img(i).mat);
            img(i).info = struct('label',label,'fov',fov,'orient',orient,'d',d,'voxsz',fov./d,'natinfo',info);
            img(i).info.voxvol = prod(img(i).info.voxsz);
            img(i).info.name = img(i).info.label;

            % Check orientation using a bone threshold
            if orientchk && opts.orient_check && check_CTlung_orientation(img(i).mat)
                writeLog(fn_log,'Permuting %s\n',tagstr{i});
                img(i).mat = permute(img(i).mat,[2,1,3]);
                img(i).info.voxsz = img(i).info.voxsz([2,1,3]);
                img(i).info.fov = img(i).info.fov([2,1,3]);
                img(i).info.d = img(i).info.d([2,1,3]);
                img(i).info.orient = img(i).info.orient([2,1,3,4],[2,1,3,4]);
            end
        end
    end

    % Check for Exp/Ins
    if any(~fnflag(:,1)) && img(1).flag && img(2).flag ...
            && check_CTlung_EIswap(img(1).mat,img(1).info.voxvol,img(2).mat,img(2).info.voxvol)
        %   ** Need to double-check in case of mislabel
        writeLog(fn_log,'Swapping INS/EXP due to lung volume\n');
        img = img([2,1]);
    end
    T.Tag(1) = {'Exp'};
    T.DataPath(1) = {img(1).sourcepath};
    T.Tag(2) = {'Ins'};
    T.DataPath(2) = {img(2).sourcepath};
    img(1).info.label = [ID,'_Exp'];
    img(2).info.label = [ID,'_Ins'];
    img(1).info.name =  [ID,'_Exp'];
    img(2).info.name =  [ID,'_Ins'];

    % Generate YACTA segmentations
    for i = 1:2
        if img(i).flag
            writeLog(fn_log,'%s : Segmentation ... ',tagstr{i});
            if fnflag(i,2)
                writeLog(fn_log,'file found: %s\n',fn{i,2});
            else
                writeLog(fn_log,'generating new ...\n');
                seg = CTlung_Segmentation(4,img(i).mat,img(i).info,img(i).info.label,procdir,fn_log);
            end

            % Check for gapped data
            if strcmp(img(i).info.natinfo.format,'DICOM') && all(isfield(img(i).info.natinfo.meta,{'SlcThk','SlicePos'})) ...
                    && (diff(img(i).info.natinfo.meta.SlicePos(1:2,3)) - img(i).info.natinfo.meta.SlcThk) > 0.0001
                oimg = img(i).mat;
                fillval = min(img(i).mat,[],'all');
                slcloc = img(i).info.natinfo.meta.SlicePos;
                % Determine gaps:
                dxyz = sqrt(sum(diff(slcloc,1).^2,2));
                d = img(i).info.d;
                if numel(unique(diff(slcloc(:,3)))) == 2
                    dnew = floor((dxyz(1)+dxyz(2))/dxyz(1));
                    dz = abs(dxyz(1)+dxyz(2))/dnew;
                    d(3) = d(3)*dnew/2;
                    ind = round(dnew/2):dnew:d(3);
                    ind = [ind;ind+1];
                else % Single-slice
                    dnew = floor(dxyz(1)/img(i).info.natinfo.meta.SlcThk);
                    dz = abs(dxyz(1))/dnew;
                    d(3) = d(3)*dnew;
                    ind = round(dnew/2):dnew:d(3);
                end
                
                disp(['Inserting image gaps. Image slices are now: ',num2str(ind(:)')])
                img(i).mat = ones(d)*fillval;
                img(i).mat(:,:,ind(:)) = oimg;
                tacq = nan(1,d(3));
                tacq(ind) = img(i).info.natinfo.meta.AcquisitionNumber;
                img(i).info.natinfo.meta.SlicePos = (0:d(3)-0.5)' * dz + slcloc(1,:);
                img(i).info.fov(3) = dz*d(3);
                img(i).info.d = d;
                img(i).info.voxsz(3) = dz;
                img(i).info.natinfo.meta.SlcThk = dz;
                img(i).info.natinfo.meta.AcquisitionNumber = tacq;
                
                % Insert gaps for the segmentation to match
                tseg = zeros(d);
                tseg(:,:,ind(:)) = seg;
                seg = tseg;
            end
            
            % Save image:
            if ~fnflag(i,1)
                saveNIFTI(fn{i,1},img(i).mat,img(i).info.label,img(i).info.fov,img(i).info.orient);
            end
            % Save segmentation
            if ~fnflag(i,2) && numel(seg)>1
                saveNIFTI(fn{i,2},seg,img(i).info.label,img(i).info.fov,img(i).info.orient);
            end
        end
    end

    dt = toc(t);
    writeLog(fn_log,'Local processing completed after %.1f minutes\n',dt/60);

    % Save results to CSV file:
    fn_res = fullfile(procdir,sprintf('%s_SourceData.csv',ID));
    writetable(T,fn_res);
catch err
    writeLog(fn_log,'ERROR in CTlung_Pipeline_local:\n%s',getReport(err,'extended','hyperlinks','off'));
    assignin('base',err);
end

