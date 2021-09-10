function [predictedLabels] = stichData_cjg_simple(net, img, classNames, inputPatchSize, outPatchSize)

% This code was pilfered from Segment3DBrainTumorUsingDeepLearningExample

% Use the overlap-tile strategy to predict the labels for each test volume.
% Each test volume is padded to make the input size a multiple of the output size of the network
% and compensates for the effects of valid convolution. The overlap-tile algorithm selects overlapping patches,
% predicts the labels for each patch by using the semanticseg function, and then recombines the patches.

tic;

disp('Processing lung volume ')

%     vol_temp = read(voldsTest); vol{id} = vol_temp(:,:,:,1); %CJG edit to 3D data (FLAIR)

% Use reflection padding for the test image.
% Avoid padding of different modalities.
volSize = size(img,(1:3));
padSizePre  = (inputPatchSize(1:3) - outPatchSize(1:3))/2;
img = padarray(img,padSizePre,'symmetric');
padSizePost = outPatchSize(1:3) - mod(volSize,outPatchSize(1:3));
img = padarray(img,padSizePost,'symmetric','post');
[heightPad,widthPad,depthPad,~] = size(img);

predictedLabels = categorical(zeros(volSize,'uint8'),(0:(numel(classNames)-1)),classNames);

% Overlap-tile strategy for segmentation of volumes.
for k = 1:outPatchSize(3):depthPad-inputPatchSize(3)+1
    for j = 1:outPatchSize(2):widthPad-inputPatchSize(2)+1
        for i = 1:outPatchSize(1):heightPad-inputPatchSize(1)+1
            patch = img( i:i+inputPatchSize(1)-1,...
                         j:j+inputPatchSize(2)-1,...
                         k:k+inputPatchSize(3)-1,:);
            patchSeg = semanticseg(patch,net);
            predictedLabels(i:i+outPatchSize(1)-1, ...
                    j:j+outPatchSize(2)-1, ...
                    k:k+outPatchSize(3)-1) = patchSeg;
            
            %CJG added to track progress
            %                 disp([i/(heightPad-inputPatchSize(1)+1) j/(widthPad-inputPatchSize(2)+1)...
            %                     k/(depthPad-inputPatchSize(3)+1)])
        end
    end
    disp(k/(depthPad-inputPatchSize(3)+1))
end

% Crop out the extra padded region.
predictedLabels = predictedLabels(1:volSize(1),1:volSize(2),1:volSize(3));

toc;