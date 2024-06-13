function [T,ver] = pipeline_vesselseg(ct,seg,info,save_path,opts_in,fn_log)
% Perform lung vessel segmentation 
% Inputs: subj_dir = directory containing all subject data
% t = vesselSeg_BH( fname_ins , fname_seg , save_path )
% t = vesselSeg_BH( ins , seg , info , save_path )

    ver = 'pipeline_vesselSeg BH-20240517';

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

    ID = info.name;
    %% Check save directory
    writeLog(fn_log,'   Saving to folder: %s\n',save_path);
    if ~isfolder(save_path)
        mkdir(save_path);
    end
    
    %% Find and load INSP CT file:
    writeLog(fn_log,'  Resampling CT image ...\n');
    voxsz_orig = info.voxsz;
    d_orig = size(ct);
    [ct,info_re] = resample_subj(ct,info,[],0.625*ones(1,3),'linear');
   
    %% Resample segmentation map:
    writeLog(fn_log,'  Resampling SEGMENTATION image ...\n');
    seg_re = resample_subj(single(seg),info,[],0.625*ones(1,3),'nearest');
    
    %% Validate size of loaded data:
    if ~all(size(ct)==size(seg_re))
        writeLog(fn_log,'Warning :: INS and SEGMENTATION images must be the same size.');
        return;
    end
    
    %% Save eroded lobe map:
    lobeval = unique(seg(seg~=0));
    fname = fullfile(save_path,[ID,'.erodedLobes.nii.gz']);
    if isfile(fname)
        writeLog(fn_log,'Loading eroded lobe map from file ...\n');
        eroded_lobes = resample_subj(readNIFTI(fname),info,info_re.d,info_re.voxsz,'nearest');
    else
        writeLog(fn_log,'Eroding lobe map ... ');
        t = tic;
        se = strel('sphere',5);
        eroded_lobes = zeros(size(seg_re));
        for i = 1:numel(lobeval)
                eroded_lobes(imerode(seg_re == lobeval(i),se)) = lobeval(i);
        end
        save_result(fname,int8(eroded_lobes),info_re,info,'nearest');
        writeLog(fn_log,'done (%s)\n',duration(0,0,toc(t)));
    end
    eroded_lobes = eroded_lobes > 0;
    
    %% Full lung volume:
    segBW = logical(seg_re);

    %% Generate enhanced vessel maps
    fname = fullfile(save_path,[ID,'.enhancedVessel.nii.gz']);
    if isfile(fname)
        writeLog(fn_log,'Loading enhanced vessel map from file ...\n');
        vessels = resample_subj(readNIFTI(fname),info,info_re.d,info_re.voxsz,'linear');
    else
        writeLog(fn_log,'Calculating enhanced vessel maps ... \n');
        t = tic;
        vessels = single(vesselSeg_subj_boxes(ct, segBW));
        vessels(isnan(vessels)) = 0;
        save_result(fname,vessels,info_re,info,'linear');
        writeLog(fn_log,'done (%s)\n\n',duration(0,0,toc(t)));
    end

    %% Binarize vessels:
    fname = fullfile(save_path,[ID,'.binVessel.nii.gz']);
    if isfile(fname)
        writeLog(fn_log,'Reading binary vessel map from file ...\n');
        bin_vessels = resample_subj(readNIFTI(fname),info,info_re.d,info_re.voxsz,'nearest');
    else
        writeLog(fn_log,'Binarizing vessel map ... ');
        t = tic;
        perLaa = nnz(ct(segBW) < -950) / nnz(segBW);
        bin_vessels = binarizeVessels(vessels,eroded_lobes,perLaa);
        save_result(fname,int8(vessels),info_re,info,'nearest');
        writeLog(fn_log,'done (%s)\n\n',duration(0,0,toc(t)));
    end
    
    %% Generate CSA maps
    fname = fullfile(save_path,[ID,'.CSA_skel.nii.gz']);
    if isfile(fname)
        writeLog(fn_log,'Reading CSA from file ...\n');
        csa_map = resample_subj(readNIFTI(fname),info,info_re.d,info_re.voxsz,'linear');
    else
        writeLog(fn_log,'Generating CSA maps ... \n');
        t = tic;
        csa_map = CSA_create_maps(bin_vessels);
        save_result(fname,csa_map,info_re,info,'linear');
        writeLog(fn_log,'done (%s)\n\n',duration(0,0,toc(t)));
    end
    
    %% Frangi Filter
    if opts.frangi
        fname = fullfile(save_path,[ID,'.enhancedVessel_frangi.nii.gz']);
        if isfile(fname)
            writeLog(fn_log,'Frangi Filtered image file found ...\n');
        else
            writeLog(fn_log,'Generating Frangi Filtered image ... \n');
            t = tic;
            opt = struct('FrangiScaleRange',[0.5 5.5],...
                         'FrangiScaleRatio',1,...
                         'BlackWhite',false);
            frangi_enhanced_vessels = FrangiFilter3D(ct.*segBW,opt) .* segBW;
            save_result(fname,frangi_enhanced_vessels,info_re,info,'linear');
            clear frangi_enhanced_vessels;
            writeLog(fn_log,'done (%s)\n\n',duration(0,0,toc(t)));
        end
    end
    
    %% Curvilinear Filter
    if opts.curvi
        fname = fullfile(save_path,[ID,'.enhancedVessel_curv.nii.gz']);
        if isfile(fname)
            writeLog(fn_log,'Curvilinear filtered image file found ...\n');
        else
            writeLog(fn_log,'Generating Curvilinear Filtered image ... \n');
            t = tic;
            curv_enhanced_vessels = vesselness3D(ct.*segBW, 0.5:1:5.5, [1,1,1], 1, true).*segBW;
            save_result(fname,curv_enhanced_vessels,info_re,info,'linear');
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
    
end

%% Resample to 0.625mm isotropic
function [I,info] = resample_subj(I,info,d_new,voxsz_new,interpm)
    if ~all(info.voxsz==voxsz_new)

        % Determine new matrix size
        if isempty(d_new)
            d_new = ceil(info.fov./voxsz_new);
        end

        dclass = class(I);
        if ~ismember(dclass,{'single','double'})
            I = single(I);
        end

        % Set up image geometry
        ext = (info.fov - info.voxsz)/2;
        F = griddedInterpolant({linspace(-ext(1),ext(1),info.d(1)),...
                                linspace(-ext(2),ext(2),info.d(2)),...
                                linspace(-ext(3),ext(3),info.d(3))},I,interpm);

        % Interpolate to new voxel locations
        ext = (d_new-1).*voxsz_new/2;
        I = F({linspace(-ext(1),ext(1),d_new(1)),...
               linspace(-ext(2),ext(2),d_new(2)),...
               linspace(-ext(3),ext(3),d_new(3))});
        
        I = cast(I,dclass);

        info.voxsz = voxsz_new;
        info.d = size(I);
    end
end

function save_result(fname,I,info,info_orig,interpm)
    if ~all(info.voxsz==info_orig.voxsz)
        I = resample_subj(I,info,info_orig.d,info_orig.voxsz,interpm);
    end
    rstr = [filesep,'([\w\.]+)\.nii'];
    if ispc
        rstr = [filesep,rstr];
    end
    label = regexp(fname,rstr,'tokens');
    if isempty(label)
        warning('Invalid file name: %s\n',fname)
    else
        saveNIFTI(fname,I,label{1},info_orig.fov,info_orig.orient);
    end
end
