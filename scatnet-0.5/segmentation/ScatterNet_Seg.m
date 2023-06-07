% Test ScatNet Segmentation
% function [atMap, laaMap, percentatMap, percentlaaMap, meanvalatMap, meanvallaaMap, miatMapLobe, milaaMapLobe, atMapLobe, laaMapLobe] = ScatNet(imgFlNm, msFlNm, flag1)
function [atMap, laaMap, percentatMap, percentlaaMap, meanvalatMap, meanvallaaMap] = ScatterNet_Seg(str, img, mask)

% Input:
%     str:            string flag for type of ScatterNet segmentation ('AT','AT_PEDS','EMPH')
%     img:            A 3D input image matrix or filename
%     mask:           Lung segmenation matrix or filename corresponding to the input image I
% Output:
%     atMap:          ScatNet map of the air trapping areas within the image
%     laaMap:         LAA map pf the air trapping areas within the image
%     percentatMap:   Percentage of ScatNet AT map
%     percentlaaMap:  Percentage of LAA AT map
%     meanvalatMap:   Mean intensity of ScatNet AT map
%     meanvallaaMap:  Mean intensity of LAA AT map

ws=25; % window size for computing features
segn=0; % number of segment. Determine automatically if set to 0
omega=.0007; % to determine the number of clusters to segment each 2D slice
var = 1.0; % to blur the image if it has a lot of noise

if ~isa(img,'char')
    I = img;
elseif contains(img,'nii')
    I = niftiread(img);
    mask = niftiread(mask);
elseif contains(img,'mhd')
    I = readMHD(img);
    mask = readMHD(mask);
end
mask = logical(mask);
np = nnz(mask);

switch str
    case 'AT'
        laaval = -856;
        fcnstr = 'FctSegAT';
    case 'AT_PEDS'
        laaval = -856;
        fcnstr = 'FctSegAT_PEDS';
    otherwise % Emph
        laaval = -950;
        fcnstr = 'FctSegEMPH';
end

roi = double(I) .* mask;
laaMap = roi <= laaval;
roi = imgaussfilt3(roi,var);
percentlaaMap = nnz(laaMap)/np;
beta = percentlaaMap;
res = feval(fcnstr,roi,ws,segn,1,omega,beta,laaval);
atMap = (res==1) & mask;

% Postprocessing
meanvallaaMap = mean(I(laaMap));
if nnz(atMap)
    atMap = logical(activecontour(roi,atMap,ceil(ws/10),'chan-vese','SmoothFactor',0,'ContractionBias',-0.1));
end
percentatMap = nnz(atMap)/np;
meanvalatMap = mean(I(atMap));


