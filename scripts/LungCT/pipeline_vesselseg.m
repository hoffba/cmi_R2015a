function [T,ver] = pipeline_vesselseg(ct,seg,info,save_path,opts_in,fn_log)
% Perform lung vessel segmentation 
% Inputs: subj_dir = directory containing all subject data
% t = vesselSeg_BH( fname_ins , fname_seg , save_path )
% t = vesselSeg_BH( ins , seg , info , save_path )

    ver = 'pipeline_vesselSeg BH-20221006';

    writeLog(fn_log,'Vessel Segmentation Version: %s\n',ver);
    
    T = [];
    tt = tic;

    %% Initialize options
    opts = struct('frangi',false,...
                  'curvi',false);
    for i = fieldnames(opts)
        if isfield(opts_in,i{1})
            opts.(i{1}) = opts_in.(i{1});
        end
    end
    
    ID = info.Description;
    %% Check save directory
    writeLog(fn_log,'   Saving to folder: %s\n',save_path);
    if ~isfolder(save_path)
        mkdir(save_path);
    end
    
    %% Find and load INSP CT file:
    writeLog(fn_log,'  Resampling CT image ...\n');
    ct = resample_subj(ct,info,fullfile(save_path,[ID,'.ct.nii']),1);
   
    %% Resample segmentation map:
    writeLog(fn_log,'  Resampling SEGMENTATION image ...\n');
    seg = resample_subj(single(seg),info,fullfile(save_path,[ID,'.lobe_segmentation.nii']),2);
    
    %% Validate size of loaded data:
    if ~all(size(ct)==size(seg))
        writeLog(fn_log,'Warning :: INS and SEGMENTATION images must be the same size.');
        return;
    end
    
    %% Update nifti info
    info.PixelDimensions = [0.625 0.625 0.625];
    info.ImageSize = size(ct);
    
    %% Save eroded lobe map:
    lobeval = unique(seg(seg~=0));
    fname = fullfile(save_path,[ID,'_erodedLobes']);
    if exist([fname,'.nii.gz'],'file')
        writeLog(fn_log,'Loading eroded lobe map from file ...\n');
        eroded_lobes = niftiread([fname,'.nii.gz']);
    else
        writeLog(fn_log,'Eroding lobe map ... ');
        t = tic;
        se = strel('sphere',5);
        eroded_lobes = zeros(size(seg));
        for i = 1:numel(lobeval)
                eroded_lobes(imerode(seg == lobeval(i),se)) = lobeval(i);
        end
        eroded_lobes = uint8(eroded_lobes);
        niftiwrite(eroded_lobes,fullfile(save_path,[ID,'_erodedLobes']),'Compressed',true);
        writeLog(fn_log,'done (%s)\n',duration(0,0,toc(t)));
    end
    eroded_lobes = eroded_lobes > 0;
    
    %% Full lung volume:
    segBW = logical(seg);

    %% Generate enhanced vessel maps
    fname = fullfile(save_path,[ID,'_enhancedVessel']);
    if exist([fname,'.nii.gz'],'file')
        writeLog(fn_log,'Loading enhanced vessel map from file ...\n');
        vessels = niftiread([fname,'.nii.gz']);
    else
        writeLog(fn_log,'Calculating enhanced vessel maps ... \n');
        t = tic;
        vessels = single(vesselSeg_subj_boxes(ct, segBW));
        vessels(isnan(vessels)) = 0;
        niftiwrite(vessels,fname,'Compressed',true);
        writeLog(fn_log,'done (%s)\n\n',duration(0,0,toc(t)));
    end

    %% Binarize vessels:
    fname = fullfile(save_path,[ID,'_binVessel']);
    if exist([fname,'.nii.gz'],'file')
        writeLog(fn_log,'Reading binary vessel map from file ...\n');
        bin_vessels = niftiread([fname,'.nii.gz']);
    else
        writeLog(fn_log,'Binarizing vessel map ... ');
        t = tic; 
        info.Datatype = 'int8';
        info.BitsPerPixel = 8;
        perLaa = nnz(ct(segBW) < -950) / nnz(segBW);
        bin_vessels = binarizeVessels(vessels,eroded_lobes,perLaa);
%         bin_vessels = int8(activecontour(vessels,imbinarize(vessels,'adaptive').*eroded_lobes,5,'Chan-Vese'));
        niftiwrite(int8(bin_vessels),fname,'Compressed',true);
        writeLog(fn_log,'done (%s)\n\n',duration(0,0,toc(t)));
    end
    
    %% Generate CSA maps
    fname = fullfile(save_path,[ID,'_CSA_skel']);
    if exist([fname,'.nii.gz'],'file')
        writeLog(fn_log,'Reading CSA from file ...\n');
        csa_map = niftiread([fname,'.nii.gz']);
    else
        writeLog(fn_log,'Generating CSA maps ... \n');
        t = tic;
        csa_map = CSA_create_maps(bin_vessels);
        niftiwrite(csa_map,fname,'Compressed',true);
        writeLog(fn_log,'done (%s)\n\n',duration(0,0,toc(t)));
    end
%     clear bin_vessels;
    
    %% Frangi Filter
    if opts.frangi
        fname = fullfile(save_path,[ID,'_enhancedVessel_frangi']);
        if exist([fname,'.nii.gz'],'file')
            writeLog(fn_log,'Frangi Filtered image file found ...\n');
        else
            writeLog(fn_log,'Generating Frangi Filtered image ... \n');
            t = tic;
            opt = struct('FrangiScaleRange',[0.5 5.5],...
                         'FrangiScaleRatio',1,...
                         'BlackWhite',false);
            frangi_enhanced_vessels = FrangiFilter3D(ct.*segBW,opt) .* segBW;
            niftiwrite(frangi_enhanced_vessels,fname,'Compressed',true);
            clear frangi_enhanced_vessels;
            writeLog(fn_log,'done (%s)\n\n',duration(0,0,toc(t)));
        end
    end
    
    %% Curvilinear Filter
    if opts.curvi
        fname = fullfile(save_path,[ID,'_enhancedVessel_curv']);
        if exist([fname,'.nii.gz'],'file')
            writeLog(fn_log,'Curvilinear filtered image file found ...\n');
        else
            writeLog(fn_log,'Generating Curvilinear Filtered image ... \n');
            t = tic;
            curv_enhanced_vessels = vesselness3D(ct.*segBW, 0.5:1:5.5, [1,1,1], 1, true).*segBW;
            niftiwrite(curv_enhanced_vessels,fname,'Compressed',true);
            clear curv_enhanced_vessels;
            writeLog(fn_log,'done (%s)\n\n',duration(0,0,toc(t)));
        end
    end
    
    %% Tabulate results and save to subject directory
    writeLog(fn_log,'Tabulating results ... ');
    t = tic;
    T = lobeLoop(seg,@(mask,binvessels,csa)vesselStats(mask,binvessels,csa),bin_vessels,csa_map);
    writetable(T,fullfile(save_path,[ID,'_vesselMetrics.csv']));
    writeLog(fn_log,'done (%s)\n\n',duration(0,0,toc(t)));
    
    writeLog(fn_log,'TOTAL vessel processing time: %s\n',duration(0,0,toc(tt)));
    
end

%% Resample to 0.625mm isotropic
function [I,info] = resample_subj(I,info,fname,tag_type)
    
%     if endsWith(fname,'.lobe_segmentation.nii.gz')
    if ~all(info.PixelDimensions==0.625)
        if tag_type == 2
            interpmethod = 'nearest';
        else
            interpmethod = 'linear';
        end
    
        %Get new size
        pixDim = info.PixelDimensions;
        pixSize = 0.625;
        
        d = size(I);
        fov = pixDim .* d;
        ext = (fov - pixDim)/2;
        F = griddedInterpolant({linspace(-ext(1),ext(1),d(1)),...
                                linspace(-ext(2),ext(2),d(2)),...
                                linspace(-ext(3),ext(3),d(3))},I,interpmethod);
        N = round(fov/pixSize);
        ext = (N-1)*pixSize/2;
        I = F({linspace(-ext(1),ext(1),N(1)),...
               linspace(-ext(2),ext(2),N(2)),...
               linspace(-ext(3),ext(3),N(3))});
        
        info.PixelDimensions = [0.625, 0.625, 0.625];
        info.ImageSize = size(I);
    end
    
    % Save resampled data to file:
    % [subj_dir,fname] = fileparts(fname);
    % fname = fullfile(subj_dir,['re_',fname,'.nii']);
    % niftiwrite(cast(I,info.Datatype),fname,info,'Compressed',true);
    
end


