function [T,ver] = pipeline_vesselseg(ct,seg,info,save_path,opts_in,fn_log)
% Perform lung vessel segmentation 
% Inputs: subj_dir = directory containing all subject data
% t = vesselSeg_BH( fname_ins , fname_seg , save_path )
% t = vesselSeg_BH( ins , seg , info , save_path )

    ver = 'pipeline_vesselSeg BH-20240617';

    T = [];
    tt = tic;

    %% Initialize options
    opts = struct('frangi',false,...
                  'curvi',false);
    if nargin>4
        for i = fieldnames(opts)
            if isfield(opts_in,i{1})
                opts.(i{1}) = opts_in.(i{1});
            end
        end
    end
    if nargin<6
        fn_log = '';
    end
    
    writeLog(fn_log,'Vessel Segmentation Version: %s\n',ver);

    ID = extractBefore(info.name,'.');

    %% Image resampling
    fn_resamp = fullfile(save_path,sprintf('%s.ins.resamp.nii.gz',ID));
    if isfile(fn_resamp)
        % Load from file
        writeLog(fn_log,'  Loading resampled image from file ...\n');
        [ct,~,fov,orient,~] = readNIFTI(fn_resamp);
        info_re = struct('fov',fov,'orient',orient,'tag','');
    else
        writeLog(fn_log,'  Resampling CT image ...\n');
        [ct,info_re] = vessel_resamp(ct,info,[],0.625*ones(1,3),'linear');
        saveNIFTI(fn_resamp,ct,'ins.resamp',info_re.fov,info_re.orient);
    end
       
    %% Resample segmentation map
    writeLog(fn_log,'  Resampling SEGMENTATION image ...\n');
    seg_re = vessel_resamp(single(seg),info,[],0.625*ones(1,3),'nearest');

    % Validate size of loaded data:
    if ~all(size(ct)==size(seg_re))
        writeLog(fn_log,'Warning :: INS and SEGMENTATION images must be the same size.');
        return;
    end

    %% Resample and erode segmentation map:    
    info_re.tag = 'vessel.erodedLobes';
    fn_erodedseg = fullfile(save_path,sprintf('%s.%s.nii.gz',ID,info_re.tag));
    if isfile(fn_erodedseg)
        writeLog(fn_log,'  Reading eroded lobes from file ...\n')
        eroded_lobes = readNIFTI(fn_erodedseg);
    else     

        % Save eroded lobe map:
        lobeval = unique(seg(seg~=0));
        writeLog(fn_log,'Eroding lobe map ... ');
        t = tic;
        se = strel('sphere',5);
        eroded_lobes = zeros(size(seg_re));
        for i = 1:numel(lobeval)
                eroded_lobes(imerode(seg_re == lobeval(i),se)) = lobeval(i);
        end
        saveNIFTI(fn_erodedseg,int8(eroded_lobes),info_re.tag,info_re.fov,info_re.orient);
        writeLog(fn_log,'done (%s)\n',duration(0,0,toc(t)));

    end
    eroded_lobes = eroded_lobes > 0;
    
    %% Full lung volume:
    segBW = logical(seg_re);

    %% Generate enhanced vessel maps
    info_re.tag = 'vessel.enhancedVessel';
    fname = fullfile(save_path,sprintf('%s.%s.nii.gz',ID,info_re.tag));
    if isfile(fname)
        writeLog(fn_log,'Loading enhanced vessel map from file ...\n');
        vessels = readNIFTI(fname);
    else
        writeLog(fn_log,'Calculating enhanced vessel maps ... \n');
        t = tic;
        vessels = single(vesselSeg_subj_boxes(ct, segBW));
        vessels(isnan(vessels)) = 0;
        saveNIFTI(fname,vessels,info_re.tag,info_re.fov,info_re.orient);
        writeLog(fn_log,'done (%s)\n\n',duration(0,0,toc(t)));
    end

    %% Binarize vessels:
    info_re.tag = 'vessel.binVessel';
    fname = fullfile(save_path,sprintf('%s.%s.nii.gz',ID,info_re.tag));
    if isfile(fname)
        writeLog(fn_log,'Reading binary vessel map from file ...\n');
        bin_vessels = readNIFTI(fname);
    else
        writeLog(fn_log,'Binarizing vessel map ... ');
        t = tic;
        perLaa = nnz(ct(segBW) < -950) / nnz(segBW);
        bin_vessels = binarizeVessels(vessels,eroded_lobes,perLaa);
        saveNIFTI(fname,int8(vessels),info_re.tag,info_re.fov,info_re.orient);
        writeLog(fn_log,'done (%s)\n\n',duration(0,0,toc(t)));
    end
    
    %% Generate CSA maps
    info_re.tag = 'vessel.CSA_skel';
    fname = fullfile(save_path,sprintf('%s.%s.nii.gz',ID,info_re.tag));
    if isfile(fname)
        writeLog(fn_log,'Reading CSA from file ...\n');
        csa_map = readNIFTI(fname);
    else
        writeLog(fn_log,'Generating CSA maps ... \n');
        t = tic;
        csa_map = CSA_create_maps(bin_vessels);
        saveNIFTI(fname,csa_map,info_re.tag,info_re.fov,info_re.orient);
        writeLog(fn_log,'done (%s)\n\n',duration(0,0,toc(t)));
    end
    
    %% Frangi Filter
    if opts.frangi
        info_re.tag = 'vessel.enhancedVessel_frangi';
        fname = fullfile(save_path,sprintf('%s.%s.nii.gz',ID,info_re.tag));
        if isfile(fname)
            writeLog(fn_log,'Frangi Filtered image file found ...\n');
        else
            writeLog(fn_log,'Generating Frangi Filtered image ... \n');
            t = tic;
            opt = struct('FrangiScaleRange',[0.5 5.5],...
                         'FrangiScaleRatio',1,...
                         'BlackWhite',false);
            frangi_enhanced_vessels = FrangiFilter3D(ct.*segBW,opt) .* segBW;
            saveNIFTI(fname,frangi_enhanced_vessels,info_re.tag,info_re.fov,info_re.orient);
            clear frangi_enhanced_vessels;
            writeLog(fn_log,'done (%s)\n\n',duration(0,0,toc(t)));
        end
    end
    
    %% Curvilinear Filter
    if opts.curvi
        info_re.tag = 'vessel.enhancedVessel_curv';
        fname = fullfile(save_path,sprintf('%s.%s.nii.gz',ID,info_re.tag));
        if isfile(fname)
            writeLog(fn_log,'Curvilinear filtered image file found ...\n');
        else
            writeLog(fn_log,'Generating Curvilinear Filtered image ... \n');
            t = tic;
            curv_enhanced_vessels = vesselness3D(ct.*segBW, 0.5:1:5.5, [1,1,1], 1, true).*segBW;
            saveNIFTI(fname,curv_enhanced_vessels,info_re.tag,info_re.fov,info_re.orient);
            clear curv_enhanced_vessels;
            writeLog(fn_log,'done (%s)\n\n',duration(0,0,toc(t)));
        end
    end
    
    %% Tabulate results and save to subject directory
    writeLog(fn_log,'Tabulating results ... ');
    t = tic;
    T = lobeLoop(seg_re,@(mask,binvessels,csa)vesselStats(mask,binvessels,csa),bin_vessels,csa_map);
    writetable(T,fullfile(save_path,[ID,'_vesselMetrics.csv']));

    

    writeLog(fn_log,'done (%s)\n\n',duration(0,0,toc(t)));
    
    writeLog(fn_log,'TOTAL vessel processing time: %s\n',duration(0,0,toc(tt)));
    

