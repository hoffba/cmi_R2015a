function T = vesselSeg_BH(varargin)
% Perform lung vessel segmentation 
% Inputs: subj_dir = directory containing all subject data
% t = vesselSeg_BH(subj_dir)
% t = vesselSeg_BH(ID,fname_ins,fname_seg,save_path)
% t = vesselSeg_BH(ID,ins,seg,info,save_path)

    T = [];
    tt = tic;

    %% Parse Inputs
    flag = true; % load file flag
    if (nargin==1) && isfolder(varargin{1})
        save_path = varargin{1};
        [~,ID] = fileparts(save_path);
        % Find files in subject directory:
        fn_ct = dir(fullfile(save_path,'*_INSP_*.ct.nii.gz'));
        if isempty(fn_ct)
            warning('No INSP file found for subject: %s.',ID);
            return
        else
            fn_ct = fullfile(save_path,fn_ct.name);
        end
        fn_seg = dir(fullfile(save_path,'*_INSP_*.lobe_segmentation.nii.gz'));
        if isempty(fn_seg)
            warning('No SEGMENTATION file found for subject: %s.',ID);
            return
        else
            fn_seg = fullfile(save_path,fn_seg.name);
        end
    elseif nargin==4
        ID = varargin{1};
        fn_ct = varargin{2};
        fn_seg = varargin{3};
        save_path = varargin{4};
    elseif nargin==5
        flag = false;
        ID = varargin{1};
        ct = varargin{2};
        seg = varargin{3};
        info = varargin{4};
        save_path = varargin{5};
    else
        error('Invalid inputs.')
    end

    %% Find and load INSP CT file:
    if flag
        if ~contains(fn_ct,'.ct.nii')
            warning('Invalid input CT file name: %s',fn_ct);
            return
        elseif exist(fn_ct,'file')
            fprintf('Reading INSP image from file: %s\n',fn_ct);
            info = niftiinfo(fn_ct);
            ct = double(niftiread(info));
            [~,fname] = fileparts(fn_ct);
            if ~startsWith(fname,'re_')
                fprintf('  Resampling CT image ...\n');
                ct = resample_subj(ct,info,fn_ct);
            end
        else
            error('Could not find file: %s',fn_ct);
        end
    elseif ~all(info.PixelDimensions==0.625)
                fprintf('  Resampling CT image ...\n');
        ct = resample_subj(ct,info,fullfile(save_path,[ID,'.ct.nii']));
    end
    
    %% Find INSP segmentation file:
    if flag
        if ~contains(fn_seg,'.lobe_segmentation.nii')
            warning('Invalid input segmentation file name: %s',fn_seg);
            return
        elseif exist(fn_seg,'file')
            fprintf('Reading SEGMENTATION image from file: %s\n',fn_seg);
            info = niftiinfo(fn_seg);
            seg = uint8(niftiread(info));
            [~,fname] = fileparts(fn_seg);
            if ~startsWith(fname,'re_')
                fprintf('  Resampling SEGMENTATION image ...\n');
                seg = resample_subj(seg,info,fn_seg);
            end
        else
            error('Could not find file: %s',fn_seg);
        end
    elseif ~all(info.PixelDimensions==0.625)
                fprintf('  Resampling SEGMENTATION image ...\n');
        seg = resample_subj(seg,info,fullfile(save_path,[ID,'.lobe_segmentation.nii']));
    end
    
    %% Update nifti info
    info.PixelDimensions = [0.625 0.625 0.625];
    info.ImageSize = size(ct);
    
    %% Save eroded lobe map:
    fprintf('Eroding lobe map ... ');
    t = tic;
    lobe = getLobeTags(seg);
    se = strel('sphere',5);
    eroded_lobes = zeros(size(seg));
    for i = 1:5
        eroded_lobes(imerode(seg == lobe(i).val,se)) = lobe(i).val;
    end
    eroded_lobes = uint8(eroded_lobes);
    save(fullfile(save_path,[ID,'_erodedLobes.mat']),'eroded_lobes');
    fprintf('done (%s)\n',duration(0,0,toc(t)));
    eroded_lobes = eroded_lobes > 0;
    
    %% Full lung volume:
    segBW = ismember(seg,[lobe.val]);

    %% Generate enhanced vessel maps
    fprintf('Calculating enhanced vessel maps ... ');
    t = tic;
    [vessels, ~] = MTHT3D( Normalize(ct.*segBW) ,...
                     0.5:1:5.5 ,...
                     12 ,... % number of orientations
                     70 ,... % beta
                     0.5 ,... % alpha
                     15 ,... % c : parameters for Vesselness
                     -1/3 );  % alfa : parameters for Neuriteness
    vessels = double(vessels) .* segBW;
    save(fullfile(save_path,[ID,'_enhancedVessel.mat']),'vessels');
    fprintf('done (%s)\n',duration(0,0,toc(t)));

    %% Binarize vessels:
    fprintf('Binarizing vessel map ... ');
    t = tic;
    info.Datatype = 'int8';
    info.BitsPerPixel = 8;
    bin_vessels = activecontour(vessels,imbinarize(vessels,'adaptive').*eroded_lobes,5,'Chan-Vese');
    niftiwrite(uint8(bin_vessels),fullfile(save_path,[ID,'_binVessel.nii']),'Compressed',true);
    fprintf('done (%s)\n',duration(0,0,toc(t)));
    
    %% Generate CSA maps
    fprintf('Generating CSA maps ... ');
    t = tic;
    csa_map = createCSA_maps(bin_vessels);
    save(fullfile(save_path,[ID,'_CSA_skel.mat']),'csa_map');
    fprintf('done (%s)\n',duration(0,0,toc(t)));
    
    %% Frangi Filter
    fprintf('Generating Frangi Filtered image ... ');
    t = tic;
    opt = struct('FrangiScaleRange',[0.5 5.5],...
                 'FrangiScaleRatio',1,...
                 'BlackWhite',false);
    frangi_enhanced_vessels = FrangiFilter3D(ct.*segBW,opt) .* segBW; 
    save(fullfile(save_path,[ID,'_enhancedVessel_frangi.mat']),'frangi_enhanced_vessels');
    fprintf('done (%s)\n',duration(0,0,toc(t)));

    %% Curvilinear Filter
    fprintf('Generating Curvilinear Filtered image ... ');
    t = tic;
    curv_enhanced_vessels = vesselness3D(ct.*segBW, 0.5:1:5.5, [1,1,1], 1, true).*segBW;
    save(fullfile(save_path,[ID,'_enhancedVessel_curv.mat']),'curv_enhanced_vessels');
    fprintf('done (%s)\n',duration(0,0,toc(t)));
    
    %% Tabulate results and save to subject directory
    fprintf('Tabulating results ... ');
    t = tic;
    T = tabulateResults(id,ct,seg,vessels,csa);
    writetable(T,fullfile(save_path,[ID,'_allMetrics.csv']));
    fprintf('done (%s)\n',duration(0,0,toc(t)));
    
    fprintf('TOTAL processing time: %s',duration(0,0,toc(tt)));
    
end

%% Resample to 0.625mm isotropic
function [I,info] = resample_subj(I,info,fname)
    
    if endsWith(fname,'.lobe_segmentation.nii.gz')
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
    [subj_dir,fname,ext] = fileparts(fname);
    fname = fullfile(subj_dir,['re_',fname,ext]);
    niftiwrite(cast(I,info.DataType),fname,info,'Compressed',true);
    
end

%% Determine lobe tags for various segmentations
function lobe = getLobeTags(seg)
    utags = unique(seg(seg>0));
    if all(ismember(1:5,utags))
        % YACTA
        lobeTag = 1:5;
        lobeName = {'LUL','LLL','RUL','RLL','RML'};
    elseif all(ismember([11,12,13,21,22],utags))
        % ImBio
        lobeTag = [11,12,13,21,22];
        lobeName = {'RUL','RLL','RML','LUL','LLL'};
    elseif all(ismember(10:10:50,utags))
        % YACTA?
        lobeTag = 10:10:50;
        lobeName = {'LUL','LLL','RUL','RLL','RML'};
    else
        error('Could not match valid segmentation labeling schema: %s',num2str(utags'));
    end
    lobe = struct('name',lobeName,'val',num2cell(lobeTag));
end

%% Tabulate statistics for each lobe
function T = tabulateResults(id,ct,seg,vessels,csa)

    lobe = getLobeTags(seg);
    nlobes = numel(lobe);
    T = table('Size',[nlobes,14],...
              'VariableTypes',{'string','string','double','double','double',...
                               'uint32','uint32','uint32',...
                               'double','double','double',...
                               'uint32','double','uint32'},...
              'VariableNames',{'ID','LOBE','VOLUME','VESSEL_VOLUME','PER_EMPH',...
                               'NUM_VESSELS','NUM_COMPONENTS','NUM_ENDPOINTS',...
                               'CSA_EXP_A','CSA_EXP_B','VESSEL_VOLUME_5DOWN',...
                               'NUM_VESSELS_5DOWN','VESSEL_VOLUME_5UP','NUM_VESSELS_5UP'});

    for i = 1:length(lobe)
        lobe_id = lobe(i).val;
        
        T.ID(i) = id;
        T.LOBE(i) = lobe(i).name;
        
        V = vessels .* (seg == lobe_id); % vessels in this lobe
        C = csa .* (seg == lobe_id); % CSA of vessels in this lobe
        
        T.VOLUME(i) = nnz(seg == lobe_id) * 0.625^3; %Convert num voxels into volume by multiplying by voxel dim
        T.VESSEL_VOLUME(i) = nnz(V) * 0.625^3;
        T.PER_EMPH(i) = nnz((ct < -950) & (seg == lobe_id)) / VOLUME*100;
       
        [T.NUM_VESSELS(i), T.NUM_COMPONENTS(i), T.NUM_ENDPOINTS(i)] = CSA_size_metrics(C);
        [T.CSA_EXP_A(i), T.CSA_EXP_B(i)] = CSA_metrics(C);
        [T.VESSEL_VOLUME_5DOWN(i), T.NUM_VESSELS_5DOWN(i)] = CSA_range_metrics(C < 5 & C > 0, C);
        [T.VESSEL_VOLUME_5UP(i), T.NUM_VESSELS_5UP(i)] = CSA_range_metrics(C > 5, C);

    end
end