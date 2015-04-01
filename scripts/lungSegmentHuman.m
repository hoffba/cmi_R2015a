% CMI script
function lungSegmentHuman(cmiObj)
%   Input: cmiObj = CMIclass object containing current settings

t = tic;

% vbody = -300; % Body threshold
% vtrach = -1024; % Trachea threshold
% vlung = -300; % Lung threshold
% se = strel('disk',2,0);

vbody = -300; % Body threshold
vtrach = -700; % Trachea threshold
vlung = -400; % Lung threshold
se = strel('disk',3,0);

cslc = cmiObj.slc(3);
cvec = cmiObj.vec;
img = cmiObj.img.mat(:,:,:,cvec);
img(img<-1024) = -1024;
d = cmiObj.img.dims(1:3);
voxvol = prod(cmiObj.img.voxsz); % in mm^3

hf = figure; colormap(gray);
hi = imagesc(img(:,:,cslc)); ha = gca; axis(ha,'off');
set(ha,'CLim',[-1024,-200]);

% First filter noise with Wiener2D:
hw = waitbar(0,'Wiener2 noise filter ...');
% for i = 1:d(3)
%     img(:,:,i) = wiener2(img(:,:,i),[3,3]);
%     img(:,:,i) = imsharpen(img(:,:,i),'Amount',2,'Radius',1);
%     waitbar(i/d(3),hw);
% end
% set(hi,'CData',img(:,:,cslc));

% Determine image ends
waitbar(0,hw,'Removing outer regions ...');
cc = bwconncomp(img<vbody);
omask = false(d);
for i = 1:cc.NumObjects
    [ii1,ii2,~] = ind2sub(d,cc.PixelIdxList{i});
    if any((ii1==1)|(ii1==d(1))|(ii2==1)|(ii2==d(2)))
        omask(cc.PixelIdxList{i}) = true;
    end
    waitbar(i/cc.NumObjects,hw);
end
img(omask) = 0;
set(hi,'CData',img(:,:,cslc));

% Find trachea
waitbar(0,hw,'Removing trachea ...');
cc = bwconncomp(img<vtrach);
ccvol = cellfun(@numel,cc.PixelIdxList)*voxvol/(10^6); % in L
tmask = false(d); tmask(cc.PixelIdxList{ccvol==max(ccvol)}) = true;
for i = 1:d(3)
    tmask(:,:,i) = imfill(imdilate(tmask(:,:,i),se),'holes');
    waitbar(i/d(3),hw);
end
img(omask|tmask) = 0;
set(hi,'CData',img(:,:,cslc));

% Find lungs
waitbar(0,hw,'Finding possible lung regions ...');
tmask = zeros(d);
for i = 1:d(3)
    
end

cc = bwconncomp(img<vlung);
ccvol = cellfun(@numel,cc.PixelIdxList)*voxvol/(10^6); % in L
ind = find(ccvol>1); nregions = length(ind);
if nregions==0
    disp('No valid regions found!');
elseif nregions>2
    disp(['Found ',num2str(nregions),' regions! = Too Many!']);
else
    se = strel('disk',3,0);
    see = strel('disk',2,0);
    mask = false([d,length(ind)]);
    for i = 1:nregions
        waitbar(0,hw,['Dilating region ',num2str(i),' of ',num2str(nregions)]);
        mask(cc.PixelIdxList{ind(i)}+(i-1)*prod(d)) = true;
        for j = 1:d(3)
            mask(:,:,j,i) = imerode(imdilate(mask(:,:,j,i),se),see);
            waitbar(j/d(3),hw);
        end
        cmiObj.img.mask.merge('Replace',mask(:,:,:,i));
        cmiObj.mask2img(['lungVOI',num2str(i)]);
    end
    cmiObj.img.mask.merge('Replace',max(mask,[],4))
    cmiObj.setSlice(round(d(3)/2));
end
delete(hw); close(hf);
disp(['Lung segmentation completed in ',num2str(toc(t)),' seconds.'])