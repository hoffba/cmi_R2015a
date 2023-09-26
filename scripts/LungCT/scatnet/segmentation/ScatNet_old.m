% Test ScatNet Segmentation
function [atMap] = ScatNet_old(imgFlNm, msFlNm, flag1)

% Input:
%     imgFlNm:      A 3D input image filename consisting of many 2D slices
%     msFlNm:       Lung segmenation mask filename corresponding to the input image I
%     flag1:        When set to 1, result of the segmentation on each 2D
%                   slice will be displayed onto the image and shown.
% Output: 
%     atMap:        label map of the air trapping areas within the image
%
%
%

ws=17; % window size for computing features
segn=0; % number of segment. Determine automatically if set to 0


I = niftiread(imgFlNm);
mask = niftiread(msFlNm);

if (length(size(I)) > 3)
    I = I(:,:,:,end);
    mask = mask(:,:,:,end);
end
   
roi = double(I).*logical(mask);
% roi = permute(roi, [2 1 3]);

mask = bwlabeln(logical(mask),26);

for i = 1:size(roi,3)

    img = roi(:,:,i);
%     figure(1), imshow(img,[]); pause(0.1); 
%     close;
    
    temp = unique(mask(:,:,i));
    if (length(temp) > 1)
    
        % Set of fixed filters
        f1=fspecial('log',[3,3],.05);
        f2=fspecial('log',[5,5],.18);
        f3=fspecial('log',[7,7],.25);

        Ig1=subImg(img,f1,f2,f3);
        Ig=cat(3,single(img),Ig1);

        res=FctSeg(Ig,ws,segn,1);
        
        % Post-Processing For Extracting the region of Air Trapping  
        
        % 1a. Extracting Left lung Clusters
        
        left = mask(:,:,i) == temp(2);
        res_left = res.*left;
        res_left = res_left + 1;
        [res_left,~] = cleanupregions(res_left,350,8);
        [res_left, ~, ~] = renumberregions(res_left);
%         res_left = res_left - 1;
        
        % 2a. Finding the unique cluster values
        
%         leftreg = 0:1:max(unique(res_left));
        
        % 3a. Computing the cluster with Minimum mean intensity value
        
        sl = struct2cell(regionprops(res_left,img,'MeanIntensity'));
%         [leftval,posleft] = min([sl{:}]);
          [posleft] = find([sl{:}] < -600);
        
        % 4a. Air Trapping Region Segmentation (Left Lung)
        
          if (isempty(posleft))
              atleft = zeros(size(img));
          else
              atleft = zeros(size(img));
              for j = 1:length(posleft)
                  atleft = atleft + (res_left == posleft(j));
              end
          end
%         if (leftval < -500)
%             atleft = res_left == leftreg(posleft+1);
%         else
%             atleft = zeros(size(img));
%         end
        
        
        % 1b. Extracting Right lung Clusters
        
        if (length(temp) > 2)
            right = mask(:,:,i) == temp(3);
            
            res_right = res.*right;
            res_right = res_right + 1;
            [res_right,~] = cleanupregions(res_right,350,8);
            [res_right, ~, ~] = renumberregions(res_right);
%             res_right = res_right - 1;
        
        % 2b. Finding the unique cluster values
        
%             rightreg = 0:1:max(unique(res_right));
            
        % 3b. Computing the cluster with Minimum mean intensity value
        
            sr = struct2cell(regionprops(res_right,img,'MeanIntensity'));
%             [rightval,posright] = min([sr{:}]);
            [posright] = find([sr{:}] < -600);  
            
        % 4b. Air Trapping Region Segmentation (Right Lung)
        
              if (isempty(posright))
              atright = zeros(size(img));
              else
                atright = zeros(size(img));
                for j = 1:length(posright)
                  atright = atright + (res_right == posright(j));
                end
              end
%             if(rightval < -500)
%                 atright = res_right == rightreg(posright+1);
%             else
%                 atright = zeros(size(img));
%             end
        else
            atright = zeros(size(img));
        end
        
        % 5. Total Air Trapping Region Segmentation
        
        atMap(:,:,i) = logical(atleft + atright);
    else
        atMap(:,:,i) = zeros(size(img));
    end
    
    % 6. show the AT regions as an overlay on the original image. 
    if (flag1 == 1)
        figure(2); imshow(imoverlay(mat2gray(img),imdilate(bwperim(atMap(:,:,i)),ones(3)),[1 0 0]),[]); pause(0.3)
    end
end
niftiwrite(im2uint8(atMap),[imgFlNm(1:end-7),'_scatnet'],'Compressed',true);
close all;


