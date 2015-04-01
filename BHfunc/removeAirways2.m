function removeAirways2(cmiObj)
% Removes airways from existing lung mask
% Assumes a filtered image is loaded (low noise from NOVA filter)

if cmiObj.img.check && cmiObj.img.mask.check
    d = cmiObj.img.dims(1:3);
    se = strel('disk',5);
    G = sobel3D(cmiObj.img.mat(:,:,:,1))<1400;
    alive = false(d);
    far = inf(d);
    
    U = 
   
    % Find initial seed point in trachea
    % * assume top of lungs is at end of slice package
    cc = bwconncomp(imclearborder(imclose(...
            cmiObj.img.mat(:,:,end,1)<-500,se)));
    np = cellfun(@length,cc.PixelIdxList);
    slcmask = false(d(1:2));
    slcmask(cc.PixelIdxList{find(np==max(np),1)}) = true;
    slcmask = imdilate(slcmask,se);
    trmask(:,:,end) = slcmask;
    ccs = bwconncomp(slcmask);
    
    
    % Remove trachea from lung mask
    trmask = cmiObj.img.mask.mat & ~trmask;
    cc = bwconncomp(trmask);
    np = cellfun(@length,cc.PixelIdxList);
    trmask(:) = false;
    trmask(cc.PixelIdxList{np==max(np)}) = true;
    cmiObj.img.mask.merge('replace',trmask);
end