% Automatic segmentation of human lung CT images
function lmask = segLungHuman_cjg_bh(output,img,dname,procdir)

%             case 2 % Clin. Insp.
%                 sm = 3;
%                 b = self.scaleB(vec);
%                 Ti = 607 + b;
%                 Tt = 125 + b;
%                 bt1 = 8; % Trachea dilation
%                 bt2 = 8; % Lung erosion
%                 bt3 = 9; % End dilation
%                 bt4 = 3; % End erosion
%                 slc = round(2*self.dims(3)/3);
% case 3 % Clin. Exp.

% CJG 20151112 added to rescale img so min = 0, this avoids having to
% do this manually using Image/Scale Image. Also changed prompt to display
% correct values of image and then scale to >0 values using "b".

if nargin<1
    output = 0; % activate prompts
end

% Initialize options:
[d1,d2,d3,d4] = size(img);
dd3 = d3*d4;
img = double(img); % convert to double



img_orig = img;
% lung_CC = 1; % default is 1 for 2nd connectivity group in lung segmentation

%% If img is dynamic (i.e. 4D)
if size(img,4)>1
    % swap phase with slices and reshape to 3D
    img = reshape(permute(img,[1 2 4 3]),d1,d2,[]);
end

% did this to remove bed from segmentation
img((end-80):end,:,:) = 0;

img_BH = img; % use this for Ben's code below my code

%% body segmentation
BW = img > -0;
se = strel('disk',10);
fprintf('Close & Fill Body:    ');
for i = 1:dd3
    fprintf('\b\b\b%u%%',round(i/dd3*100));
    BW(:,:,i) = imerode(imfill(imdilate(BW(:,:,i),se)),se);
end
fprintf('\n');
CC = bwconncomp(BW);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idx] = max(numPixels);
BW = false(d1,d2,dd3);
BW(CC.PixelIdxList{idx}) = true;

%% median filter [9 9]
filt_img = img;
fprintf('Medfilt2 No Parallel:    ')
for i = 1:dd3
    fprintf('\b\b\b%u%%',round(i/dd3*100));
    filt_img(:,:,i) = medfilt2(img(:,:,i).*BW(:,:,i),[9 9]);
end
fprintf('\n');

%% Lung and airway segmentation
BW = (filt_img <= -250) & BW;
fprintf('ImClearBorder:    ');
for i = 1:size(BW,3)
    fprintf('\b\b\b%u%%',round(i/dd3*100));
    BW(:,:,i) = imclearborder(BW(:,:,i));
end
fprintf('\n');
CC = bwconncomp(BW);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idx] = max(numPixels);
BW = false(d1,d2,dd3);
BW(CC.PixelIdxList{idx}) = true;

se = strel('disk',5);
if d4 > 1
    sl_ph = d4*(d3-1)+1; % first slice of last phase
    mask_temp = imerode(imfill(BW(:,:,sl_ph)),strel('disk',6));
    
    CC = bwconncomp(mask_temp);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [~,idx] = min(numPixels);
    
    mask = false(size(img));
    
    mask_temp2 = false(size(mask_temp));
    mask_temp2(CC.PixelIdxList{idx}) = 1;
    
    mask(:,:,sl_ph) = logical(imdilate(mask_temp2,strel('disk',9)));
    
    MedTrach = 74; % default for trachea threshold 74=-950; 124=-900
else
    trach_size = sum(reshape(BW(:,:,:),[],size(BW,3)),1);
    
    for CCi = size(BW,3):-1:1
        CC = bwconncomp(BW(:,:,CCi));
        
        if CC.NumObjects > 1
            break;
            
        end
    end
    trach_size(1,1:CCi) = 0;
    
    trach_idx = find(sum(trach_size,1)>200&sum(trach_size,1)<2000);
    if isempty(trach_idx)
        trach_idx = size(BW,3);
    end
    
    mask = false(size(filt_img));
    for tIdx = 1: length(trach_idx)
        mask(:,:,trach_idx(tIdx)) = logical(imdilate(BW(:,:,trach_idx(tIdx)),se)); % identify trachea for Ben code
    end
    
%     img_trach = img(:,:,trach_idx); BW_trach = BW(:,:,trach_idx);
%     MedTrach = prctile(img_trach(logical(BW_trach))+1024,[75]);

    MedTrach = 74; % default for trachea threshold
end

%% Start human segmentation
img = img_BH; % use original unfiltered img

if min(reshape(img,[],1))<0
    img=img+1024;
    b=-1024;
else
    b=0; %min
end

sm = 3;
Ti = 774 + b; % 774=-250, 874=-150

% Tt = 124 + b;
Tt = MedTrach + b; % calculated the median HU in trachea above

bt1 = 8; % Trachea dilation
bt2 = 8; % Lung erosion
bt3 = 9; % End dilation
bt4 = 3; % End erosion
d = size(img);
slc = round(2*d(3)/3);
dopt = ~isempty(slc);

answer = {num2str(sm), num2str(Ti), num2str(Tt), num2str(bt1), num2str(bt2),...
          num2str(bt3), num2str(bt4)};
if output == 0
    prompt = {'Smoothing','Lung Threshold','Trachea Threshold',...
        'Trachea Dilation','Lung Erosion','End Dilation','End Erosion'};
    answer = inputdlg(prompt,'Input Parameters',1,answer);
end

if ~isempty(answer)
    val = str2double(answer{1});
    if ~isnan(val)
        sm = val;
    end
    val = str2double(answer{2});
    if ~isnan(val)
        Ti = val-b;
    end
    val = str2double(answer{3});
    if ~isnan(val)
        Tt = val-b;
    end
    val = str2double(answer{4});
    if ~isnan(val)
        bt1 = val;
    end
    val = str2double(answer{5});
    if ~isnan(val)
        bt2 = val;
    end
    val = str2double(answer{6});
    if ~isnan(val)
        bt3 = val;
    end
    val = str2double(answer{7});
    if ~isnan(val)
        bt4 = val;
    end
    
    disp('.')
    disp('.')
    disp('Starting Lung Segmentation ...')
    tic
    
    if sm
        img = filtGaussSep(img,sm);
    end
    % Make voxels outside FOV have air HU value
    img(img<b) = b;
    
    % Find image ends
    tmask = false(d);
    tmask([1,end],:,:) = true;
    tmask(:,[1,end],:) = true;
    cind = find(tmask);
    vind = find(mask);
    
    % Create initial filtered masks for lung and trachea
    if dopt
        th = tic;
    end
    % Trachea mask
    for j = 1:size(img,3)
        lmask(:,:,j) = filtGaussSep(double(img(:,:,j)>=Tt),2) <= 0.1; % erode to get rid of small regions
        lmask(:,:,j) = filtGaussSep(lmask(:,:,j),bt1) > 0.1; % dilate to fill trachea
    end
    
%     lmask = BW_body2; % used this to replace above code.
    
    cc = bwconncomp(lmask);
    lmask = false(d);
    for i = 1:cc.NumObjects
        if any(ismember(vind,cc.PixelIdxList{i}))
            lmask(cc.PixelIdxList{i}) = true;
        end
    end
    
    se = strel('disk',5);
    for j = 1:size(lmask,3)
        lmask(:,:,j) = filtGaussSep(double(lmask(:,:,j)),bt1) > 0.01; % dilate
    end
    
    %     for j = 1:size(lmask,3)
    %         lmask(:,:,j) = imdilate(lmask(:,:,j), se);
    %     end
    
    if dopt
        hf_segm=figure('Position',[40 40 600 800]);
%         hf_segm = figure(1000); % CJG 20210512 did this to prevent trying to close Reg when open
        subplot(3,2,1),
        imshow(lmask(:,:,slc));
        title('Trachea')
        disp(['  Filterd trachea mask (' num2str(toc(th)) ' seconds)...'])
        pause(0.1)
        th = tic;
    end
    % Lung mask - w/o trachea
    
    for j = 1:size(lmask,3)
        lmask(:,:,j) = filtGaussSep(1-double((img(:,:,j)<Ti) & ~lmask(:,:,j)),bt2) < 0.1; % erode
    end
    
    %     for j = 1:size(lmask,3)
    %         lmask(:,:,j) = imerode(1-double(lmask(:,:,j)), se);
    %     end
    
    h100 = waitbar(0,'ImClearBorder');
    for i = 1:size(BW,3)
        waitbar(i/size(img,3),h100,[str,'Clear Borders: ',num2str(round(100*i/size(img,3))),' % completed']);
        lmask(:,:,i) = imclearborder(lmask(:,:,i));
    end
    close(h100)
    if dopt
        figure(hf_segm),subplot(3,2,2),
        imshow(lmask(:,:,slc));
        title('Lungs')
        disp(['  Filtered lung mask (' num2str(toc(th)) ' seconds)...'])
        pause(0.1)
        th = tic;
    end
    
    % Find connected regions
    cc = bwconncomp(lmask);
    if dopt
        disp(['  Found ' num2str(cc.NumObjects) ' regions (' num2str(toc(th)) ' seconds)...'])
    end
    
    % Discard regions that are too small and largest (outside of body)
    n = numel(lmask);
    nvox = cellfun(@numel,cc.PixelIdxList);
    ind = cellfun(@(x)any(ismember(x,cind)),cc.PixelIdxList);
    
    clear idx
    %     idx = find(~ind & ((nvox/n)>0.002)); % CJG this threshold desides to keep or disgard volume
    
    a = sort(nvox,'descend');
    idx(1) = find((nvox == a(1)));
    idx(2) = find((nvox == a(2)));
    
    disp(['    Using ' num2str(length(idx)) ' regions ...'])
    %     lmask2 = false(d); % temporary mask
    lmask = zeros(d);  % final label matrix
    %     lmask_label = zeros(d); % CJG 20160304 generate label for individual masks
    
    se = strel('disk',10);
    for i = 1:length(idx)
        disp(['    Region ' num2str(i) ': ' num2str(nvox(idx(i))) ' voxels'])
        if dopt
            th2 = tic;
        end
        lmask2 = false(d);
        lmask2(cc.PixelIdxList{idx(i)}) = true;
        
        %         lmask2 = filtGaussSep(double(lmask2),bt3) > 0.1; % dilate
        
        for j = 1:size(lmask2,3)
            lmask2(:,:,j) = imdilate(lmask2(:,:,j), se);
        end
        
        if dopt
            figure(hf_segm),subplot(3,2,2*i+1),
            imshow(lmask2(:,:,slc));
            title('Dilated')
            disp(['      Dilated (' num2str(toc(th2)) ' seconds)...'])
            pause(0.1)
            th2 = tic;
        end
        cc2 = bwconncomp(~lmask2);
        nvox2 = cellfun(@numel,cc2.PixelIdxList);
        lmask2 = true(d);
        lmask2(cc2.PixelIdxList{nvox2==max(nvox2)}) = false;
        if dopt
            disp(['      Found outside (' num2str(toc(th2)) ' seconds)...'])
            th2 = tic;
        end
        
        %         lmask2 = filtGaussSep(1-double(lmask2),bt4) < 0.1; % erode
        
        for j = 1:size(lmask2,3)
            lmask2(:,:,j) = imerode(lmask2(:,:,j), se);
        end
        if dopt
            figure(hf_segm),subplot(3,2,2*i+2),
            imshow(lmask2(:,:,slc));
            title('Eroded')
            disp(['      Eroded (' num2str(toc(th2)) ' seconds)...'])
            pause(0.1)
        end
        
        lmask(lmask2) = i;
        %         lmask = (lmask | lmask2);
        %         lmask_temp(:,:,:,i)=lmask2.*i; % CJG 20160304 used for generating label
        %         lmask_label=sum(lmask_temp,4); % CJG 20160304 sum lung masks.
    end
    
    lmask = int8(lungSubDiv(lmask));
    disp(['Ended after ' num2str(toc) ' seconds'])
    
    %%
%     [PathName a z] = fileparts(info.Filename);
%     [x a z] = fileparts(a);
%     
%     fNiftiname = fullfile(PathName,[a,'.nii']);
%     
%     info.Filename = fNiftiname;
%     info.Datatype = 'int8';
%     info.BitsPerPixel = 8;
    
    %% If img is dynamic (i.e. 4D)
    if size(img_orig,4)>1
        % reshape to 4D [rows columns phase slice] and permute phase and slice
        lmask = permute(reshape(lmask,size(img_orig,1), size(img_orig,2), size(img_orig,4), size(img_orig,3)),[1 2 4 3]);
    end
    
%     ['Output  :',a,'_Label_Results.jpg']
    if nargin==4 && isfolder(procdir)
        saveas(hf_segm,fullfile(procdir,[dname,'_Results.jpg'])); 
        close(hf_segm)
%     niftiwrite(int8(lmask),fNiftiname,info,'Compressed',true);
%     cd(PathName);
    end
end