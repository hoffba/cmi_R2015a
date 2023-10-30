function brainSeg_MatlabEx

modelFolder = fileparts(mfilename("fullpath"));
modelFile = fullfile(modelFolder,"trainedBrainSynthSegNetwork.h5");

dataDir = "R:\CGalban_Lab\LabMembers\BenHoff\tempDATA\DIPG\brainSegData";

imFile = fullfile(dataDir,"anat.nii.gz");
metaData = niftiinfo(imFile);
vol = niftiread(metaData);

labelDirs = fullfile(dataDir,"groundTruth");
[classNames,labelIDs] = getBrainCANDISegmentationLabels;

resample = false;
cropSize = 192;
[volProc,cropIdx,imSize] = preProcessBrainCANDIData(vol,metaData,cropSize,resample);
inputSize = size(volProc);

volDL = dlarray(volProc,"SSSCB");

lgraph = importKerasLayers(modelFile,ImportWeights=true,ImageInputSize=inputSize);

placeholderLayers = findPlaceholderLayers(lgraph);
sf = softmaxLayer;
for i = 1:numel(placeholderLayers)
    lgraph = replaceLayer(lgraph,placeholderLayers(i).Name,sf);
end
net = dlnetwork(lgraph);
layerGraph(net)

predictIm = predict(net,volDL);
flipVal = true;
if flipVal
    flippedData = fliplr(volProc);  
    flippedData = flip(flippedData,2);
    flippedData = flip(flippedData,1);
    flippedData = dlarray(flippedData,"SSSCB");
    flipPredictIm = predict(net,flippedData);
else
    flipPredictIm = [];  
end

predictedSegMaps = postProcessBrainCANDIData(predictIm,flipPredictIm,imSize,cropIdx,metaData,classNames,labelIDs);
sliceIdx = 80;
testSlice = rescale(volProc(:,:,sliceIdx));
predSegMap = predictedSegMaps(:,:,sliceIdx);
B = labeloverlay(testSlice,predSegMap,"IncludedLabels",2:32);
figure
montage({testSlice,B})

