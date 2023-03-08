% Test ScatNet Segmentation
% function [atMap, laaMap, percentatMap, percentlaaMap, meanvalatMap, meanvallaaMap, miatMapLobe, milaaMapLobe, atMapLobe, laaMapLobe] = ScatNet(imgFlNm, msFlNm, flag1)
function [atMap, laaMap, percentatMap, percentlaaMap, meanvalatMap, meanvallaaMap] = ScatterNetAT(imgFlNm, msFlNm, flag1)

% Input:
%     imgFlNm:      A 3D input image filename consisting of many 2D slices
%     msFlNm:       Lung segmenation mask filename corresponding to the input image I
%     flag1:        When set to 1, result of the segmentation on each 2D
%                   slice will be displayed onto the image and shown.
%                   (Just for visualization.)
% Output: 
%     atMap:          ScatNet map of the air trapping areas within the image
%     laaMap:         LAA map pf the air trapping areas within the image
%     percentatMap:   Percentage of ScatNet AT map
%     percentlaaMap:  Percentage of LAA AT map
%     meanvalatMap:   Mean intensity of ScatNet AT map
%     meanvallaaMap:  Mean intensity of LAA AT map
%     miatMapLobe:    Mean intensity of ScatNet AT map for each lobe
%     milaaMapLobe:   Mean intensity of LAA AT map for each lobe
%     atMapLobe:      Percentage of ScatNet AT map for each lobe
%     laaMapLobe:     Percentage of LAA AT map for each lobe
%
%
%
ws=25; % window size for computing features
segn=0; % number of segment. Determine automatically if set to 0
omega=.0007; % to determine the number of clusters to segment each 2D slice
var = 1.0; % to blur the image if it has a lot of noise
laaval = -856;

% NIFTIREAD for Nifti file formats
I = niftiread(imgFlNm);
mask = niftiread(msFlNm);

% %MHDREAD for MHD file formats
% [I,label,fov,info] = readMHD(imgFlNm);
% [mask,label2,fov2,info2] = readMHD(msFlNm);

 if (length(size(I)) > 3)
    I = I(:,:,:,end);
    mask = mask(:,:,:,end);
end
   
roi = double(I).*logical(mask);
% roi = permute(roi, [2 1 3]);
mask2 = mask;
mask = ((mask == 1) | (mask == 2));
% mask = bwlabeln(logical(mask),26);

% Calculate the -856HU Threshold and the %AT due to the LAA Method
laaMap = roi <= laaval;

percentlaaMap = length(find(laaMap))/length(find(mask));
beta = percentlaaMap;
meanvallaaMap = mean(I(laaMap));

res = FctSegAT(imgaussfilt3(roi,var),ws,segn,1,omega,beta,laaval);

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

            % 3a. Air Trapping Region Segmentation (Left Lung)

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

            % 3b. Air Trapping Region Segmentation (Right Lung)

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

            % 4. Total Air Trapping Region Segmentation

            atMap(:,:,i) = logical(atleft + atright);
        else
            atMap(:,:,i) = zeros(size(img));
        end
    else
        atMap(:,:,i) = logical((res(:,:,i) == 1) .* logical(mask2(:, :, i)));
        
    end
    % 5. show the AT regions as an overlay on the original image. 
    if (flag1 == 1)
        figure(2); imshow(imoverlay(mat2gray(img'),imdilate(bwperim(atMap(:,:,i)'),ones(3)),[1 0 0]),[]); pause(0.0000000001)
%         temp1 = imoverlay(mat2gray(img'),imdilate(bwperim(atMap(:,:,i)'),ones(2)),[1 0 0]);
%         temp2 = imoverlay(mat2gray(img'),imdilate(bwperim(laaMap(:,:,i)'),ones(2)),[0 1 0]);
%         imwrite(temp1,'atMapOV.tif','tif','WriteMode','append');
%         imwrite(temp2,'laaMapOV.tif','tif','WriteMode','append');
    end
end
% Postprocessing
atMap = activecontour(imgaussfilt3(roi,var),atMap,ceil(ws/10),'chan-vese','SmoothFactor',0,'ContractionBias',-0.1);
% avgatMap = mean(I(atMap));
% atMap = gaussianBlur(roi,var)<= avgatMap;
percentatMap = length(find(atMap))/length(find(mask2));
meanvalatMap = mean(I(atMap));

% % Calculate the mean Intensity values for each lobe
% val1 = mean(I(logical(atMap)));
% val2 = mean(I(logical(laaMap)));
% miatMapLobe = [];
% milaaMapLobe = [];
% r1 = unique(mask2(:));
% for i = 1:length(r1)
%     temp = mean(I(logical((atMap.*logical(mask2 == r1(i))))));
%     temp2 = mean(I(logical((laaMap.*logical(mask2 == r1(i))))));
%     miatMapLobe = [miatMapLobe, temp];
%     milaaMapLobe = [milaaMapLobe, temp2];
% end
% 
% %Calculate percent AT Map values for each lobe
% atMapLobe = [];
% laaMapLobe = [];
% for i = 1:length(r1)
% %   t1 = length(find(logical((atMap.*logical(mask2 == r1(i))))))/length(find(mask));
%     t2 = length(find(logical((atMap.*logical(mask2 == r1(i))))))/length(find(logical(mask2 == r1(i))));
% %   t3 = length(find(logical((laaMap.*logical(mask2 == r1(i))))))/length(find(mask));
%     t4 = length(find(logical((laaMap.*logical(mask2 == r1(i))))))/length(find(logical(mask2 == r1(i))));
%     atMapLobe = [atMapLobe, t2];
%     laaMapLobe = [laaMapLobe, t4];
% end

% niftiwrite(im2uint8(atMap),[imgFlNm(1:end-7),'_scatternetat'],'Compressed',true);
% niftiwrite(im2uint8(atMap),[imgFlNm(1:end-4),'_ScatterNetAT'],'Compressed',true);
close all;

