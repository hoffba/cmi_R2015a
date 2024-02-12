function [seg,vols] = brainSeg_SynthSeg(varargin)
% [seg,vols] = brainSeg_SynthSeg(nifti_filename)
%   - Performs segmentation on image in input Nifti file
%   - Writes result to file, with *.SynthSeg.* tag
% [seg,vols] = brainSeg_SynthSeg(img,info)
%   - Performs segmentation on input image matrix
%   - Requires image geometry paramters in the info structure (in Nifti format)
%       * Fields: .Transform.T ; .PixelDimensions
% [classNames,labelIDs] = brainSeg_SynthSeg();
%   - Returns segmentation labeling scheme
%
% **SynthSeg: Segmentation of brain MRI scans of any contrast and resolution without retraining** \
% B. Billot, D.N. Greve, O. Puonti, A. Thielscher, K. Van Leemput, B. Fischl, A.V. Dalca, J.E. Iglesias \
% Medical Image Analysis (2023) \
% [ [article](https://www.sciencedirect.com/science/article/pii/S1361841523000506) | [arxiv](https://arxiv.org/abs/2107.09559) | [bibtex](bibtex.bib) ]

seg = []; vols = [];
fname = '';
write_flag = false;
qc_chk = false;

%% Handle inputs
if nargin==0
    % Return labeling scheme if no inputs
    [classNames,labelIDs] = getBrainCANDISegmentationLabels();
    seg = classNames;
    vols = labelIDs;
    return;
elseif ischar(varargin{1})
    fname = varargin{1};
    if nargin==2
        qc_chk = varargin{2};
    end
elseif (nargin>1) && isnumeric(varargin{1}) && isstruct(varargin{2})
    img = varargin{1};
    info = varargin{2};
    if nargin==3
        qc_chk = varargin{3};
    end
end

t = tic;

%% Load image from file
if fname
    % Load image from file
    if isfolder(fname)
        [img,label,fov,orient,~] = cmi_load(1,[],fname);
        d = size(img,1:3);
        info = init_niftiinfo(label,fov./d,'single',d,orient);
        img = single(img);
    elseif isfile(fname)
        info = niftiinfo(fname);
        img = niftiread(info);
        write_flag = true;
    else
        return;
    end
end

%% Find image geometry
dims = size(img);
aff = info.Transform.T';
voxsz = info.PixelDimensions;

%% Preprocess image
    disp('... Pre-processing segmentation map')
% Resample image
    [X1,aff] = resampleData(img,voxsz,aff);% Align volume to ref
% Align to reference geometry
    [X1,~] = AlignToRef(X1, aff, eye(4));
% pad if needed
    numLevels = 5;
    N = 2^numLevels;
    dims_ref = size(X1);
    sizNew = N*ceil(dims_ref/N);
    dataNew = zeros(sizNew,class(X1));
    dataNew(1:dims_ref(1),1:dims_ref(2),1:dims_ref(3)) = X1;
    X1 = dataNew; clear dataNew
% Normalize
    inMax = prctile(X1(:),99.5);
    inMin = prctile(X1(:),0.5);
    X1 = rescale(X1,0,1,InputMax=inMax,InputMin=inMin);

%% Find and load UNet model file:
modelFolder = fileparts(mfilename("fullpath"));
modelFile = fullfile(modelFolder,"trainedBrainSynthSegNetwork.h5");
if ~isfile(modelFile)
    warning('UNet model file not found: %s',modelFile);
    return;
end
warning('off','all')
lgraph = importKerasLayers(modelFile,ImportWeights=true,ImageInputSize=size(X1));
% placeholderLayers = findPlaceholderLayers(lgraph);
sf = softmaxLayer;
lgraph = replaceLayer(lgraph,"unet_prediction",sf);
net = dlnetwork(lgraph);
warning('on','all')

%% Run model prediction(s)
% Primary prediction
    disp('... Predicting segmentation in original orientation')
    X2 = dlarray(X1,"SSSCB");
    warning('off','all')
    predictIm = predict(net,X2);
    warning('on','all')
% Run flipped prediction
    flip_flag = false;
    if flip_flag
        disp('... Predicting segmentation in flipped orientation')
        fX2 = dlarray(flip(X1,1),"SSSCB");
        warning('off','all')
        flipPredictIm = predict(net,fX2);
        warning('on','all')
    else
        flipPredictIm = [];
    end

%% Postprocess the segmentation
    disp('... Post-processing segmentation map')
% Extract and smooth primary seg data
    predictIm = squeeze(extractdata(predictIm));
    predictIm = predictIm(1:dims_ref(1), 1:dims_ref(2), 1:dims_ref(3), :);
    sigma = 0.5;
    for n = 1:size(predictIm,4)
        predictIm(:,:,:,n) = imgaussfilt3(predictIm(:,:,:,n),sigma,FilterSize=3,Padding=0);   
    end
% Extract, smooth, and correct labels for flipped seg data
    if ~isempty(flipPredictIm)
        flipPredictIm = squeeze(extractdata(flipPredictIm));
        flipPredictIm = flipPredictIm(1:dims_ref(1), 1:dims_ref(2), 1:dims_ref(3), :);
        % Apply 3D Gaussian smoothing on fliped data
        for n = 1:size(flipPredictIm,4)
            flipPredictIm(:,:,:,n) = imgaussfilt3(flipPredictIm(:,:,:,n),sigma,FilterSize=3,Padding=0);
        end
        % Get posteriors and segmentation
        flipPredictIm = flip(flipPredictIm,1);
        lrIndices = [2  3  4  5  6  7  8  9  10 11 15 16 17 18; ...
                     19 20 21 22 23 24 25 26 27 28 29 30 31 32];
        rlIndices = reshape(flip(lrIndices)',1,[]);
        lrIndices = reshape(lrIndices',1,[]);
        flipPredictIm(:,:,:,lrIndices) = flipPredictIm(:,:,:,rlIndices);
        predictIm = 0.5 * (predictIm + flipPredictIm);
    end
% Use largest connected component
    th = 0.25;
    temp = predictIm(:,:,:,2:end);
    postInMask = sum(temp,4)>th;
    largeComp = findLargestComponent(postInMask);
    S = repmat(largeComp,1,1,1,size(temp,4));
    temp(~S)=0;
    predictIm(:,:,:,2:end) = temp;
% Make posteriors to zero outside the largest connected component of each topological class
    postInMask = predictIm > th;
    topology_classes = [0  4  4  5  5  6  6  7  8  9 10  1  2  3  5  11 ...
                        12 13 14 14 15 15 16 16 17 18 19 20 15 21 22 23];
    topology_classesUnique = unique(topology_classes);
    for topology_class = topology_classesUnique(2):topology_classesUnique(end)
        [~,tmp_topology_indices] = find(topology_classes==topology_class);
        tmp_mask = postInMask(:,:,:,tmp_topology_indices);
        tmp_mask =  findLargestComponent(tmp_mask);
        predictIm(:,:,:,tmp_topology_indices) = predictIm(:,:,:,tmp_topology_indices).*tmp_mask;
    end
% Renormalize posteriors and get hard segmentation
    predictIm = predictIm./sum(predictIm,4);
    [~,seg] = max(predictIm,[],4);
    [classNames,labelIDs] = getBrainCANDISegmentationLabels;
    seg = uint8(labelIDs(seg)); 
% align back to native space
    [seg, ~] = AlignToRef(seg, eye(4), aff);
% Interpolate back to original
    seg = resampleData(seg,voxsz,dims);

%% Calculate segmentation region volumes
    if nargout == 2
        vols = table('Size',[numel(classNames),2],'VariableTypes',{'cellstr','double'},...
            'VariableNames',{'ROI','Volume_mm3'});
        vols.Volume_mm3 = squeeze(sum(sum(sum(predictIm,1),2),3));
        vols.ROI = classNames';
    end

%% Show a few slices
    if logical(qc_chk)
        labelSlice3D(img,seg,voxsz);
    end

%% Write segmenation to file
    if write_flag
        gzflag = endsWith(fname,'.gz');
        [fpath,fname] = fileparts(fname);
        fname = fullfile(fpath,[extractBefore(fname,'.nii'),'.SynthSeg.nii']);
        disp('... Writing segmentation to file');
        info.Datatype = 'uint8';
        info.BitsPerPixel = 8;
        niftiwrite(seg, fname, info, Compressed=gzflag);
    end

fprintf('Segmentation complete (%s)\n',string(datetime(0,0,0,0,0,toc(t),'Format','mm:ss')));

function [B,aff] = resampleData(A,voxsz,x)
% Interpolates:
% 1) Image to isotropic (1mm)^3 using linear interpolation
% 2) Segmentation back to original geometry using nearest neighbor
d = size(A);
class_str = class(A);
imflag = all(size(x)==[4,4]);
if imflag
    % Interpolate image
    aff = x;
    interpmethod = 'linear';
    voxsz_in = voxsz;
    voxsz_out = ones(1,3);
    % Update affine transform
    aff(1:3,1:3)= aff(1:3,1:3)./voxsz_in;
    aff(1:3,end) = aff(1:3,end)-(aff(1:3,1:3)*0.5*(voxsz_in'-1));
    % Determine new matrix size
    N = ceil(voxsz.*d);
else
    % Intepolate segmentation
    interpmethod = 'nearest';
    voxsz_in = ones(1,3);
    voxsz_out = voxsz;
    N = x;
end
ext = (d-1).*voxsz_in/2;
F = griddedInterpolant({linspace(-ext(1),ext(1),d(1)),...
                        linspace(-ext(2),ext(2),d(2)),...
                        linspace(-ext(3),ext(3),d(3))},single(A),interpmethod);
ext = (N-1).*voxsz_out/2;
B = F({linspace(-ext(1),ext(1),N(1)),...
       linspace(-ext(2),ext(2),N(2)),...
       linspace(-ext(3),ext(3),N(3))});
B = cast(B,class_str);


function largeCompFinal =  findLargestComponent(postInMask)
largeCompFinal = zeros(size(postInMask));
for ii = 1:size(postInMask,4)
    tmp =  postInMask(:,:,:,ii);
    if any(tmp,'all')
        % get_largest_connected_component
        CC = bwconncomp(postInMask(:,:,:,ii),6);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [~,idx] = max(numPixels);
        largeComp = false(size(tmp));
        largeComp(CC.PixelIdxList{idx}) = true;
        largeCompFinal(:,:,:,ii) = largeComp;
    end
end


function [classNames,labelIDs] = getBrainCANDISegmentationLabels()
% Define segmentation class names and numeric label IDs
classNames = ["background","leftCerebralWhiteMatter","leftCerebralCortex",...
    "leftLateralVentricle","leftInferiorLateralVentricle","leftCerebellumWhiteMatter",...
    "leftCerebellumCortex","leftThalamus","leftCaudate","leftPutamen","leftPallidum",...
    "thirdVentricle","fourthVentricle","brainStem","leftHippocampus","leftAmygdala",...
    "leftAccumbensArea","leftVentralDC","rightCerebralWhiteMatter","rightCerebralCortex",...
    "rightLateralVentricle","rightInferiorLateralVentricle","rightCerebellumWhiteMatter",...
    "rightCerebellumCortex","rightThalamus","rightCaudate","rightPutamen","rightPallidum",...
    "rightHippocampus","rightAmygdala","rightAccumbensArea","rightVentralDC"];

labelIDs = [0 2 3 4 5 7 8 10 11 12 13 14 15 16 17 18 26, ...
    28 41 42 43 44 46 47 49 50 51 52 53 54 58 60];

