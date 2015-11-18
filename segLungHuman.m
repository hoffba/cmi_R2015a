% Automatic segmentation of human lung CT images
function lmask = segLungHuman(img,mask)

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
if min(img(mask))<0
    img=img+1024;
    b=-1024;
else
    b=0; %min
end
    
sm = 3;
Ti = 607 + b;
Tt = 125 + b;
bt1 = 8; % Trachea dilation
bt2 = 8; % Lung erosion
bt3 = 9; % End dilation
bt4 = 3; % End erosion
slc = round(2*size(img,3)/3);
dopt = ~isempty(slc);

prompt = {'Smoothing','Lung Threshold','Trachea Threshold',...
          'Trachea Dilation','Lung Erosion','End Dilation','End Erosion'};
def = {num2str(sm),num2str(Ti),num2str(Tt),num2str(bt1),num2str(bt2),...
       num2str(bt3),num2str(bt4)};
answer = inputdlg(prompt,'Input Parameters',1,def);
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
    tmask = false(size(img));
    tmask([1,end],:,:) = true;
    tmask(:,[1,end],:) = true;
    cind = find(tmask);
    vind = find(mask);

    % Create initial filtered masks for lung and trachea
    if dopt
        th = tic;
    end
    % Trachea mask
    lmask = filtGaussSep(double(img>=Tt),2) <= 0.1; % erode to get rid of small regions
    lmask = filtGaussSep(lmask,bt1) > 0.1; % dilate to fill trachea

    cc = bwconncomp(lmask);
    lmask = false(size(lmask));
    for i = 1:cc.NumObjects
        if any(ismember(vind,cc.PixelIdxList{i}))
            lmask(cc.PixelIdxList{i}) = true;
        end
    end

    lmask = filtGaussSep(double(lmask),bt1) > 0.01; % dilate
    if dopt
        hf=figure('Position',[40 40 600 800]);
        subplot(3,2,1),
        imshow(lmask(:,:,slc));
        title('Trachea')
        disp(['  Filterd trachea mask (' num2str(toc(th)) ' seconds)...'])
        pause(0.1)
        th = tic;
    end
    % Lung mask - w/o trachea
    lmask = filtGaussSep(1-double((img<Ti) & ~lmask),bt2) < 0.1; % erode
    if dopt
        figure(hf),subplot(3,2,2),
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
    idx = find(~ind & ((nvox/n)>0.005));
    disp(['    Using ' num2str(length(idx)) ' regions ...'])
    lmask2 = false(size(lmask)); lmask = lmask2;
    for i = 1:length(idx)
        disp(['    Region ' num2str(i) ': ' num2str(nvox(idx(i))) ' voxels'])
        if dopt
            th2 = tic;
        end
        lmask2 = false(size(lmask));
        lmask2(cc.PixelIdxList{idx(i)}) = true;
        lmask2 = filtGaussSep(double(lmask2),bt3) > 0.1; % dilate
        if dopt
            figure(hf),subplot(3,2,2*i+1),
            imshow(lmask2(:,:,slc));
            title('Dilated')
            disp(['      Dilated (' num2str(toc(th2)) ' seconds)...'])
            pause(0.1)
            th2 = tic;
        end
        cc2 = bwconncomp(~lmask2);
        nvox2 = cellfun(@numel,cc2.PixelIdxList);
        lmask2 = true(size(lmask));
        lmask2(cc2.PixelIdxList{nvox2==max(nvox2)}) = false;
        if dopt
            disp(['      Found outside (' num2str(toc(th2)) ' seconds)...'])
            th2 = tic;
        end
        lmask2 = filtGaussSep(1-double(lmask2),bt4) < 0.1; % erode
        if dopt
            figure(hf),subplot(3,2,2*i+2),
            imshow(lmask2(:,:,slc));
            title('Eroded')
            disp(['      Eroded (' num2str(toc(th2)) ' seconds)...'])
            pause(0.1)
        end
        lmask = (lmask | lmask2);
    end
    disp(['Ended after ' num2str(toc) ' seconds'])
end