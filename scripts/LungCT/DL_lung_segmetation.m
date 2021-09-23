function finalLabel = DL_lung_segmetation(x)

if nargin<1
    [FileName,PathName] = uigetfile('*.nii.gz','Select nii File for Conversion of mhd 2 nii.');
    if FileName
        img = niftiread(fullfile(PathName, FileName));
        info = niftiinfo(fullfile(PathName, FileName));
        dims = info.ImageSize;
    end
elseif ischar(x) && exist(x,'file')
    img = permute(niftiread(x), [2 1 3]);
    info = niftiinfo(x);
    dims = info.ImageSize;
elseif isnumeric(x) && ndims(x)==3
    img = x;
    dims = size(x);
else
    error('Invalid input.');
end
clear x;

%% Pre-processing:
img(img<-1024) = -1024; % get rid of values <-1024
img = medfilt3(img, [3 3 3]);
% Downsize image matrix for processing
procdims = 128.*ones(1,3);
img = single(imresize3(img,procdims,'method','nearest'));

%% Generate body mask
% Threshold softtissue
finalLabel = img>-250;
finalLabel(:,:,end) = imfill(finalLabel(:,:,end),'holes'); % remove trachea (disconnect from boundary)
finalLabel(:,:,1) = imfill(finalLabel(:,:,1), 'holes'); % this helps capture all low density objects in body
finalLabel = imclose(finalLabel,strel('cube',3));

% need to fill holes on sides in case CT is off center and lung
% hits boundary
finalLabel = permute(finalLabel, [1 3 2]);
finalLabel(:,:,1) = imfill(finalLabel(:,:,1), 'holes');
finalLabel(:,:,end) = imfill(finalLabel(:,:,end), 'holes');
finalLabel = permute(finalLabel, [1 3 2]);

finalLabel = imfill(finalLabel,'holes'); % fill entire image

% Keep only largest connected component
CC = bwconncomp(finalLabel);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idx] = max(numPixels);
finalLabel = false(procdims);
finalLabel(CC.PixelIdxList{idx}) = true;

%% Normalize the image
img = (single(img.*finalLabel) + 1024) / (max(img(finalLabel)) + 1024); % cjg data must be in range [0 1]

%% Load CNN: CJG got validation accuracy to 98.8% Still get some colon and trachea. Use post-processing to remove.
net = load(fullfile(fileparts(which('DL_lung_segmetation')),...
    'trained3DUNetValid-01-Sep-2021-10-52-44-Epoch-60.mat'));
net = net.net;

%% Stich Data
classNames = ["background","Lungs"];
inputPatchSize = [116 116 116];
outPatchSize = [28 28 28 2];
lungLabel = stichData_cjg_simple(net, img, classNames, inputPatchSize, outPatchSize);

%% upsample back to original dimensions
lungLabel = double(imresize3(single(lungLabel),dims,'nearest'));
lungLabel = lungLabel-1;

%% VOI post-processing
% erode image
se = strel('cube',5);
finalLabel = imerode(lungLabel,se);

% find lungs: sometimes colon is captured
finalLabel = bwlabeln(finalLabel);
A = unique(finalLabel(:,:,round(dims(3)/2)));
A(A==0)=[];
if isempty(A)
    warning('No lung VOI found.');
    finalLabel = lungLabel;
else
    fprintf('Found %u lung regions.\n',numel(A));
    finalLabel(~ismember(finalLabel,A)) = 0;
    for i = 1:numel(A)
        finalLabel(finalLabel==A(i)) = i;
    end
    % dilate image
    finalLabel = imdilate(finalLabel,se);
end


