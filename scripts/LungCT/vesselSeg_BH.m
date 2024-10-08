function [T,ver] = vesselSeg_BH(varargin)
% Perform lung vessel segmentation 
% Inputs: subj_dir = directory containing all subject data
% t = vesselSeg_BH( fname_ins , fname_seg , save_path )
% t = vesselSeg_BH( ins , seg , info , save_path )

    ver = 'vesselSeg_BH SR-20220623';

    T = [];
    tt = tic;

    fn_ct = '';
    fn_seg = '';
    ct = [];
    seg = [];
    info = [];
    
    %% Parse Inputs
    flag = true; % load file flag
    if nargin==3
        fn_ct = varargin{1};
        fn_seg = varargin{2};
        save_path = varargin{3};
        
        [~,ID] = fileparts(fn_ct);
        if startsWith(ID,'re_')
            ID(1:3) = [];
        end
        if contains(ID,'_')
            ID = extractBefore(ID,'_');
        elseif contains(ID,'.')
            ID = extractBefore(ID,'.');
        end
    elseif nargin==4
        ct = varargin{1};
        seg = varargin{2};
        info = varargin{3};
        save_path = varargin{4};
        
        ID = info.Description;
    else
        error('Invalid inputs.')
    end

    %% Check save directory
    fprintf('   Saving to folder: %s\n',save_path);
    if ~isfolder(save_path)
        mkdir(save_path);
    end
    
    %% Find and load INSP CT file:
    tag_type = 1;
    if isempty(ct) % filename input
        if ~contains(fn_ct,'.nii')
            warning('Invalid input CT file name: %s',fn_ct);
            return
        elseif exist(fn_ct,'file')
            fprintf('Reading INSP image from file: %s\n',fn_ct);
            info = niftiinfo(fn_ct);
            ct = double(niftiread(info));
            [~,fname] = fileparts(fn_ct);
            if ~startsWith(fname,'re_')
                fprintf('  Resampling CT image ...\n');
                ct = resample_subj(ct,info,fn_ct,tag_type);
            end
        else
            error('Could not find file: %s',fn_ct);
        end
    elseif ~all(info.PixelDimensions==0.625)
                fprintf('  Resampling CT image ...\n');
        ct = resample_subj(ct,info,fullfile(save_path,[ID,'.ct.nii']),tag_type);
    end
   
    %% Find INSP segmentation file:
    tag_type = 2;
    if isempty(seg)
        if ~contains(fn_seg,'.nii')
            warning('Invalid input segmentation file name: %s',fn_seg);
            return
        elseif exist(fn_seg,'file')
            fprintf('Reading SEGMENTATION image from file: %s\n',fn_seg);
            info = niftiinfo(fn_seg);
            seg = uint8(niftiread(info));
            [~,fname] = fileparts(fn_seg);
            if ~startsWith(fname,'re_')
                fprintf('  Resampling SEGMENTATION image ...\n');
                seg = resample_subj(seg,info,fn_seg,tag_type);
            end
        else
            error('Could not find file: %s',fn_seg);
        end
    elseif ~all(info.PixelDimensions==0.625)
                fprintf('  Resampling SEGMENTATION image ...\n');
        seg = resample_subj(seg,info,fullfile(save_path,[ID,'.lobe_segmentation.nii']),tag_type);
    end
    
    %% Validate size of loaded data:
    if ~all(size(ct)==size(seg))
        warning('Warning :: INS and SEGMENTATION images must be the same size.');
        return;
    end
    
    %% Update nifti info
    info.PixelDimensions = [0.625 0.625 0.625];
    info.ImageSize = size(ct);
    
    %% Save eroded lobe map:
    lobeval = unique(seg(seg~=0));
    fname = fullfile(save_path,[ID,'_erodedLobes']);
    if flag && exist([fname,'.nii.gz'],'file')
        fprintf('Loading eroded lobe map from file ...\n');
        eroded_lobes = niftiread([fname,'.nii.gz']);
    else
        fprintf('Eroding lobe map ... ');
        t = tic;
        se = strel('sphere',5);
        eroded_lobes = zeros(size(seg));
        for i = 1:numel(lobeval)
                eroded_lobes(imerode(seg == lobeval(i),se)) = lobeval(i);
        end
        eroded_lobes = uint8(eroded_lobes);
        niftiwrite(eroded_lobes,fullfile(save_path,[ID,'_erodedLobes']),'Compressed',true);
        fprintf('done (%s)\n',duration(0,0,toc(t)));
    end
    eroded_lobes = eroded_lobes > 0;
    
    %% Full lung volume:
    segBW = logical(seg);

    %% Generate enhanced vessel maps
    fname = fullfile(save_path,[ID,'_enhancedVessel']);
    if flag && exist([fname,'.nii.gz'],'file')
        fprintf('Loading enhanced vessel map from file ...\n');
        vessels = niftiread([fname,'.nii.gz']);
    else
        fprintf('Calculating enhanced vessel maps ... \n');
        t = tic;
        vessels = single(vesselSeg_subj_boxes(ct, segBW));
        vessels(isnan(vessels)) = 0;
        niftiwrite(vessels,fname,'Compressed',true);
        fprintf('done (%s)\n\n',duration(0,0,toc(t)));
    end

    %% Binarize vessels:
    fname = fullfile(save_path,[ID,'_binVessel']);
    if flag && exist([fname,'.nii.gz'],'file')
        fprintf('Reading binary vessel map from file ...\n');
        bin_vessels = niftiread([fname,'.nii.gz']);
    else
        fprintf('Binarizing vessel map ... ');
        t = tic; 
        info.Datatype = 'int8';
        info.BitsPerPixel = 8;
        perLaa = nnz(ct(segBW) < -950) / nnz(segBW);
        bin_vessels = binarizeVessels(vessels,eroded_lobes,perLaa);
%         bin_vessels = int8(activecontour(vessels,imbinarize(vessels,'adaptive').*eroded_lobes,5,'Chan-Vese'));
        niftiwrite(int8(bin_vessels),fname,'Compressed',true);
        fprintf('done (%s)\n\n',duration(0,0,toc(t)));
    end
    
    %% Generate CSA maps
    fname = fullfile(save_path,[ID,'_CSA_skel']);
    if flag && exist([fname,'.nii.gz'],'file')
        fprintf('Reading CSA from file ...\n');
        csa_map = niftiread([fname,'.nii.gz']);
    else
        fprintf('Generating CSA maps ... \n');
        t = tic;
        csa_map = CSA_create_maps(bin_vessels);
        niftiwrite(csa_map,fname,'Compressed',true);
        fprintf('done (%s)\n\n',duration(0,0,toc(t)));
    end
%     clear bin_vessels;
    
    %% Frangi Filter
%     fname = fullfile(save_path,[ID,'_enhancedVessel_frangi']);
%     if flag && exist([fname,'.nii.gz'],'file')
%         fprintf('Frangi Filtered image file found ...\n');
%     else
%         fprintf('Generating Frangi Filtered image ... \n');
%         t = tic;
%         opt = struct('FrangiScaleRange',[0.5 5.5],...
%                      'FrangiScaleRatio',1,...
%                      'BlackWhite',false);
%         frangi_enhanced_vessels = FrangiFilter3D(ct.*segBW,opt) .* segBW;
%         niftiwrite(frangi_enhanced_vessels,fname,'Compressed',true);
%         clear frangi_enhanced_vessels;
%         fprintf('done (%s)\n\n',duration(0,0,toc(t)));
%     end

    %% Curvilinear Filter
%     fname = fullfile(save_path,[ID,'_enhancedVessel_curv']);
%     if flag && exist([fname,'.nii.gz'],'file')
%         fprintf('Curvilinear filtered image file found ...\n');
%     else
%         fprintf('Generating Curvilinear Filtered image ... \n');
%         t = tic;
%         curv_enhanced_vessels = vesselness3D(ct.*segBW, 0.5:1:5.5, [1,1,1], 1, true).*segBW;
%         niftiwrite(curv_enhanced_vessels,fname,'Compressed',true);
%         clear curv_enhanced_vessels;
%         fprintf('done (%s)\n\n',duration(0,0,toc(t)));
%     end
    
    %% Tabulate results and save to subject directory
    fprintf('Tabulating results ... ');
    t = tic;
    T = lobeLoop(seg,@(mask,ins,binvessels,csa)vesselStats(mask,ins,binvessels,csa),ct,bin_vessels,csa_map);
    writetable(T,fullfile(save_path,[ID,'_vesselMetrics.csv']));
    fprintf('done (%s)\n\n',duration(0,0,toc(t)));
    
%     T = lobeTable2struct(T);
%     ID = {ID};
%     T = addvars(T,ID,'Before',1);
    
    fprintf('TOTAL processing time: %s',duration(0,0,toc(tt)));
    
end

%% Resample to 0.625mm isotropic
function [I,info] = resample_subj(I,info,fname,tag_type)
    
%     if endsWith(fname,'.lobe_segmentation.nii.gz')
    if tag_type == 2
        interpmethod = 'nearest';
    else
        interpmethod = 'linear';
    end

    %Get new size
    pixDim = info.PixelDimensions;
    pixSize = 0.625;
    
    [sx,sy,sz] = size(I);
    xq = length(0:pixSize/pixDim(1):sx);
    yq = length(0:pixSize/pixDim(2):sy);
    zq = sz; %was length(0:sz) so original batch has one extra z layer
    
    I = imresize3(I, [xq yq zq], 'method', interpmethod);
    
    info.PixelDimensions = [0.625, 0.625, 0.625];
    info.ImageSize = size(I);
    
    % Save resampled data to file:
    [subj_dir,fname] = fileparts(fname);
    fname = fullfile(subj_dir,['re_',fname,'.nii']);
    niftiwrite(cast(I,info.Datatype),fname,info,'Compressed',true);
    
end


