function seg = brainSeg(imFile)
% From Matlab example

flag_flip = true;

%% Set up model
modelFolder = fileparts(mfilename("fullpat"));
modelFile = fullfile(modelFolder,"trainedBrainSynthSegNetwork.h5");
if ~isfile(modelFile)
    trainedBrainCANDINetwork_url = "https://www.mathworks.com/supportfiles/image/data/trainedBrainSynthSegNetwork.h5";
    downloadTrainedNetwork(trainedBrainCANDINetwork_url,modelFolder);
end

flag_p = false;
if nargin
    [procdir,id] = fileparts(imFile);
    id = extractBefore(id,'.');
    flag_p = true;
else
    % Matlab example
    procdir = "R:\CGalban_Lab\LabMembers\BenHoff\tempDATA\DIPG";
    dataDir = fullfile(procdir,"brainSegData");
    imFile = fullfile(dataDir,"anat.nii.gz");
    labelDirs = fullfile(dataDir,"groundTruth");
    if ~isfile(imFile)
        zipFile = matlab.internal.examples.downloadSupportFile("image","data/brainSegData.zip");
        unzip(zipFile,procdir)
    end
end

% Read image and info
info = niftiinfo(imFile);
vol = single(niftiread(info));
voxsz_orig = info.PixelDimensions;

% Adjust image orientation based on model
% Model need RAS coordinate system
if flag_p
    p = [3,1,2];
    vol = flip(permute(vol,p),1);
    voxsz_orig = voxsz_orig(p);
else
    p = [1,3,2]; f = [3];
    vol = permute(vol,p);
    for i = 1:numel(f)
        vol = flip(vol,f(i));
    end
    voxsz_orig = voxsz_orig(p);
end

% Identify brain classes in model
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

%% Preprocess data
% resample image volume to isotropic voxel size (1mm^3)
voxsz_re = ones(1,3);
matsz_re = 192*ones(1,3);
matsz_orig = size(vol,1:3);
volProc = resampleImg(vol,voxsz_orig,voxsz_re,matsz_re,"linear");
volProc = int16(volProc);
% Normalize image values
inMax = prctile(volProc(:),99.5);
inMin = prctile(volProc(:),0.5);
volProc = rescale(volProc,0,1,InputMax=inMax,InputMin=inMin);
% Rotate image to reference space

% convert image to DL array for processing
volDL = dlarray(volProc,"SSSCB");

%% Set up DL model
lgraph = importKerasLayers(modelFile,ImportWeights=true,ImageInputSize=matsz_re);
sf = softmaxLayer;
lgraph = replaceLayer(lgraph,"unet_prediction",sf);
net = dlnetwork(lgraph);
layerGraph(net)

%% Run DL prediction
predictIm = predict(net,volDL);
if flag_flip
    flippedData = fliplr(volProc);  
    flippedData = flip(flippedData,2);
    flippedData = flip(flippedData,1);
    flippedData = dlarray(flippedData,"SSSCB");
    flipPredictIm = predict(net,flippedData);
end

%% Postprocess segmentation map
% Apply 3D Gaussian smoothing
predictIm = squeeze(extractdata(predictIm));
sigma = 0.5;
for n = 1:size(predictIm,4)
    predictIm(:,:,:,n) = imgaussfilt3(predictIm(:,:,:,n),sigma,FilterSize=3,Padding=0);   
end
% Process flipped result
if flag_flip
    flipPredictIm = squeeze(extractdata(flipPredictIm));
    % Apply 3D Gaussian smoothing on fliped data
    for n = 1:size(flipPredictIm,4)
        flipPredictIm(:,:,:,n) = imgaussfilt3(flipPredictIm(:,:,:,n),sigma,FilterSize=3,Padding=0);
    end
    % Get posteriors and segmentation
    flipPredictIm = flip(flipPredictIm,1);
    flipPredictIm = flip(flipPredictIm,2);
    flipPredictIm = fliplr(flipPredictIm);
    lrIndices = [2 3 4 5 6 7 8 9 10 11 15 16 17 18; 19 20 21 22 ...
        23 24 25 26 27 28 29 30 31 32];
    rlIndices = flip(lrIndices);
    flipPredictIm(:,:,:,reshape(lrIndices',1,[])) = flipPredictIm(:,:,:,reshape(rlIndices',1,[]));
    predictIm = (predictIm + flipPredictIm)/2;
end
% Use largest connected component
th = 0.25;
temp = predictIm(:,:,:,2:end);
postInMask = sum(temp,4)>th;
largeCompFinal = findLargestComponent(postInMask);
S = repmat(largeCompFinal,1,1,1,size(temp,4));
temp(~S)=0;
predictIm(:,:,:,2:end) = temp;
% Make posteriors to zero outside the largest connected component 
% of each topological class
postInMask = predictIm > th;
topology_classes = [0  4  4  5  5  6  6  7  8  9 10  1  2  3  5 11 ...
    12 13 14 14 15 15 16 16 17 18 19 20 15 21 22 23];
topology_classesUnique = unique(topology_classes);
for topology_class = topology_classesUnique(2):topology_classesUnique(end)
    [~,tmp_topology_indices] = find(topology_classes==topology_class);
    tmp_mask = postInMask(:,:,:,tmp_topology_indices);
    if any(tmp_mask,'all')
        tmp_mask =  findLargestComponent(tmp_mask);
        predictIm(:,:,:,tmp_topology_indices) = predictIm(:,:,:,tmp_topology_indices).*tmp_mask;
    end
end
% Renormalize posteriors and get hard segmentation
predictIm = predictIm./sum(predictIm,4);
[~,predictedSegMaps] = max(predictIm,[],4);


% Display middle slice
predictedSegMaps = categorical(predictedSegMaps,labelIDs,classNames);
displayOverlay(volProc,predictedSegMaps)


% Re-map segmentation results to original geometry
predictedSegMaps = resampleImg(predictedSegMaps,voxsz_re,voxsz_orig,matsz_orig,"nearest");
% Convert segmentation result to categorical type
predictedSegMaps = categorical(predictedSegMaps,labelIDs,classNames);

%% Save result to .nii
info.Datatype = 'uint8';
info.BitsPerPixel = 8;
niftiwrite(uint8(predictedSegMaps),fullfile(procdir,[id,'.MatlabSeg']),'Compressed',true);


function I = resampleImg(I,voxsz_in,voxsz_out,matsz_out,interp_method)
    % Input matrix size
    matsz_in = size(I,1:3);
    % Set up interpolant
    ext = voxsz_in.*(matsz_in-1)/2;
    F = griddedInterpolant({linspace(-ext(1),ext(1),matsz_in(1)),...
                            linspace(-ext(2),ext(2),matsz_in(2)),...
                            linspace(-ext(3),ext(3),matsz_in(3))},...
                            I, interp_method, "none");
    % Interpolate to desired voxel locations
    ext = voxsz_out.*(matsz_out-1)/2;
    I = F({linspace(-ext(1),ext(1),matsz_out(1)),...
           linspace(-ext(2),ext(2),matsz_out(2)),...
           linspace(-ext(3),ext(3),matsz_out(3))});

    I(isnan(I)) = 0;

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

function displayOverlay(volProc,predictedSegMaps)
    sliceIdx = round(size(volProc,3)/2);
    testSlice = rescale(volProc(:,:,sliceIdx));
    predSegMap = predictedSegMaps(:,:,sliceIdx);
    B = labeloverlay(testSlice,predSegMap,"IncludedLabels",2:32);
    figure
    montage({testSlice,B})

