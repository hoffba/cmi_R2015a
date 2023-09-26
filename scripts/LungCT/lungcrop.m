function timg = lungcrop(img,voxvol)

slc_pad = 10;

% Binarize image
timg = medfilt2_3D(img);
timg = timg < -150;

% Remove anything connected to the image edges in 2D
for i = 1:size(timg,3)
    slc = timg(:,:,i);
    cc = bwconncomp(slc);
    for j = 1:cc.NumObjects
        [ii,jj] = ind2sub(cc.ImageSize,cc.PixelIdxList{j});
        if any(ismember([1,cc.ImageSize(1)],ii)) || any(ismember([1,cc.ImageSize(2)],jj))
            slc(cc.PixelIdxList{j}) = false;
        end
    end
    timg(:,:,i) = slc;
end

% Largest connected region is likely the lungs
cc = bwconncomp(timg);
N = cellfun(@numel,cc.PixelIdxList);
timg(:) = false;
timg(cc.PixelIdxList{N==max(N)}) = true;

% Find top and bottom of lungs
V = squeeze(sum(timg,[1,2]))*voxvol;
mxi = find(V==max(V),1);
bottomi = max(find(V,1)-slc_pad,1);
topi = min(find(V(mxi:end)<1000,1)+mxi-1+slc_pad,cc.ImageSize(3)); % threshold @ 1000ul

% Crop original image
timg = img(:,:,bottomi:topi);