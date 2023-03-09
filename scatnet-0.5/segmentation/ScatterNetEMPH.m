% Test ScatNet Segmentation
% function [emphMap, laaMap, percentemphMap, percentlaaMap, meanvalemphMap, meanvallaaMap, miemphMapLobe, milaaMapLobe, emphMapLobe, laaMapLobe] = ScatterNetEMPH(imgFlNm, msFlNm, flag1)
function [emphMap, laaMap, percentemphMap, percentlaaMap, meanvalemphMap, meanvallaaMap] = ScatterNetEMPH(imgFlNm, msFlNm, flag1)

% Input:
%     imgFlNm:      A 3D input image filename consisting of many 2D slices
%     msFlNm:       Lung segmenation mask filename corresponding to the input image I
%     flag1:        When set to 1, result of the segmentation on each 2D
%                   slice will be displayed onto the image and shown.
% Output: 
%     emphMap:          ScatNet map of the air trapping areas within the image
%     laaMap:           LAA map of the emphysematous areas within the image
%     percentemphMap:   Percentage of ScatNet EMPH map
%     percentlaaMap:    Percentage of LAA EMPH map
%     meanvalatMap:     Mean intensity of ScatNet EMPH map
%     meanvallaaMap:    Mean intensity of LAA EMPH map
%     miemphMapLobe:    Mean intensity of ScatNet EMPH map for each lobe
%     milaaMapLobe:     Mean intensity of LAA EMPH map for each lobe
%     emphMapLobe:      Percentage of ScatNet EMPH map for each lobe
%     laaMapLobe:       Percentage of LAA EMPH map for each lobe
%
%
%
ws=25; % window size for computing features
segn=0; % number of segment. Determine automatically if set to 0
omega=.0007; % to determine the number of clusters to segment each 2D slice
var = 1.0; % to blur the image if it has a lot of noise
laaval = -950;

% % NIFTIREAD for Nifti file formats
% I = niftiread(imgFlNm);
% mask = niftiread(msFlNm);
% 
% % %MHDREAD for MHD file formats
% % [I,label,fov,info] = readMHD(imgFlNm);
% % [mask,label2,fov2,info2] = readMHD(msFlNm);

if ~isa(imgFlNm,'char')
    I = imgFlNm;
    mask = msFlNm;
    
elseif contains(imgFlNm,'nii')
    % NIFTIREAD for Nifti file formats
    I = niftiread(imgFlNm);
    mask = niftiread(msFlNm);
    
elseif contains(imgFlNm,'mhd')
    % %MHDREAD for MHD file formats
    I = readMHD(imgFlNm);
    mask = readMHD(msFlNm);
    
end

 if (length(size(I)) > 3)
    I = I(:,:,:,end);
    mask = mask(:,:,:,end);
 end
   
d = size(I);
 
roi = double(I).*logical(mask);
% roi = permute(roi, [2 1 3]);
mask2 = mask;
mask = ((mask == 1) | (mask == 2));
% mask = bwlabeln(logical(mask),26);

% Calculate the -950HU Threshold and the %EMPH due to the LAA Method
laaMap = roi <= laaval;

percentlaaMap = length(find(laaMap))/length(find(mask));
beta = percentlaaMap;
meanvallaaMap = mean(I(laaMap));

% res = FctSegEmph(gaussianBlur(roi,var),ws,segn,1,omega,beta);
res = FctSegEMPH(imgaussfilt3(roi,var),ws,segn,1,omega,beta,laaval);

emphMap = false(d);
for i = 1:size(roi,3)

    img = roi(:,:,i);
    if (segn ~= 0)
    %     figure(1), imshow(img,[]); pause(0.1); 
    %     close;

        temp = unique(mask(:,:,i));
        if (length(temp) > 1)

            % Set of fixed filters
            f1=fspecial('log',[3,3],.05);
            f2=fspecial('log',[5,5],.18);
            f3=fspecial('log',[7,7],.25);
            f4=fspecial('log',[9,9],.35);
            f5=gabor_fn(1.5,pi/2);
            f6=gabor_fn(1.5,0);
            f7=gabor_fn(1.5,pi/4);
            f8=gabor_fn(1.5,-pi/4);
            f9=gabor_fn(2.5,pi/2);
            f10=gabor_fn(2.5,0);
            f11=gabor_fn(2.5,pi/4);
            f12=gabor_fn(2.5,-pi/4);

            Ig1=subImg(img,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12);
            Ig=cat(3,single(img),Ig1);

            res=double(FctSeg(Ig,ws,segn,1,omega)).*logical(mask(:,:,i));
    %         figure(1); imshow(res',[]);
            % Post-Processing For Extracting the region of Air Trapping  

            %Constraining using Mean intensity of the region
            meanInt = -600;

            % 1a. Extracting Left lung Clusters

            left = mask(:,:,i) == temp(2);
            res_left = res.*left;
            res_left = res_left + 1;
            [res_left,~] = cleanupregions(res_left,150,8);
            [res_left, ~, ~] = renumberregions(res_left);
    %         figure(1); subplot 121; imshow(res_left',[]);
            % 2a. Finding the unique cluster values

            sl = struct2cell(regionprops(res_left,img,'MeanIntensity'));
    %         [leftval,posleft] = min([sl{:}]);
              [posleft] = find([sl{:}] < meanInt);

            % 3a. Emphysema Region Segmentation (Left Lung)

              if (isempty(posleft))
                  atleft = zeros(size(img));
              else
                  atleft = zeros(size(img));
                  for j = 1:length(posleft)
                      atleft = atleft + (res_left == posleft(j));
                  end
              end

            % 1b. Extracting Right lung Clusters

            if (length(temp) > 2)
                right = mask(:,:,i) == temp(3);

                res_right = res.*right;
                res_right = res_right + 1;
                [res_right,~] = cleanupregions(res_right,150,8);
                [res_right, ~, ~] = renumberregions(res_right);
    %             res_right = res_right - 1;
    %                 subplot 122; imshow(res_right',[]);
            % 2b. Finding the unique cluster values  

                sr = struct2cell(regionprops(res_right,img,'MeanIntensity'));
    %             [rightval,posright] = min([sr{:}]);
                [posright] = find([sr{:}] < meanInt);  

            % 3b. Emphysema Region Segmentation (Right Lung)

                  if (isempty(posright))
                    atright = zeros(size(img));
                  else
                    atright = zeros(size(img));
                    for j = 1:length(posright)
                      atright = atright + (res_right == posright(j));
                    end
                  end
            else
                atright = zeros(size(img));
            end

            % 4. Total Emphysema Region Segmentation

            emphMap(:,:,i) = logical(atleft + atright);
        end
    else
        emphMap(:,:,i) = logical((res(:,:,i) == 1) .* logical(mask2(:, :, i)));
    end
    % 5. show the AT regions as an overlay on the original image. 
    if (flag1 == 1)
        figure(2); imshow(imoverlay(mat2gray(img'),imdilate(bwperim(emphMap(:,:,i)'),ones(3)),[1 0 0]),[]); pause(0.0000000001)
%         temp1 = imoverlay(mat2gray(img'),imdilate(bwperim(emphMap(:,:,i)'),ones(2)),[1 0 0]);
%         temp2 = imoverlay(mat2gray(img'),imdilate(bwperim(laaMap(:,:,i)'),ones(2)),[0 1 0]);
%         imwrite(temp1,'emphMapOV.tif','tif','WriteMode','append');
%         imwrite(temp2,'laaMapOV.tif','tif','WriteMode','append');
    end
end
% Postprocessing
if nnz(emphMap)
    emphMap = activecontour(imgaussfilt3(roi,var),emphMap,ceil(ws/10),'chan-vese','SmoothFactor',0,'ContractionBias',-0.1);
    % avgemphMap = mean(I(emphMap));
    % emphMap = gaussianBlur(roi,var)<= avgemphMap;
end
percentemphMap = nnz(emphMap)/nnz(mask2);
meanvalemphMap = mean(I(emphMap));
    

% % Calculate the mean Intensity values for each lobe
% val1 = mean(I(logical(emphMap)));
% val2 = mean(I(logical(laaMap)));
% miemphMapLobe = [];
% milaaMapLobe = [];
% r1 = unique(mask2(:));
% for i = 1:length(r1)
%     temp = mean(I(logical((emphMap.*logical(mask2 == r1(i))))));
%     temp2 = mean(I(logical((laaMap.*logical(mask2 == r1(i))))));
%     miemphMapLobe = [miemphMapLobe, temp];
%     milaaMapLobe = [milaaMapLobe, temp2];
% end
% 
% %Calculate percent EMPH Map values for each lobe
% emphMapLobe = [];
% laaMapLobe = [];
% for i = 1:length(r1)
% %   t1 = length(find(logical((emphMap.*logical(mask2 == r1(i))))))/length(find(mask));
%     t2 = length(find(logical((emphMap.*logical(mask2 == r1(i))))))/length(find(logical(mask2 == r1(i))));
% %   t3 = length(find(logical((laaMap.*logical(mask2 == r1(i))))))/length(find(mask));
%     t4 = length(find(logical((laaMap.*logical(mask2 == r1(i))))))/length(find(logical(mask2 == r1(i))));
%     emphMapLobe = [emphMapLobe, t2];
%     laaMapLobe = [laaMapLobe, t4];
% end

% niftiwrite(im2uint8(emphMap),[imgFlNm(1:end-7),'_scatternetemph'],'Compressed',true);
% niftiwrite(im2uint8(emphMap),[imgFlNm(1:end-4),'_ScatterNetEMPH'],'Compressed',true);
close all;


