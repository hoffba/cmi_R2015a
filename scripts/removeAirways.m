function removeAirways(cmiObj)
% Removes airways from existing lung mask
% Assumes a filtered image is loaded (low noise from NOVA filter)

if cmiObj.img.check && cmiObj.img.mask.check
    % Set growing parameters:
    T = -700; % Image threshold
    A = 1.5; % 2D region size threshold
    se = strel('disk',5);
    Tmask = min(cmiObj.img.mat(:,:,:,1),[],4) < T;
    G = sobel3D(cmiObj.img.mat(:,:,:,1)) < 1400;
    
    % Initialize 3D matrices:
    d = cmiObj.img.dims(1:3);
    trmask = false(d);
    
    % Find initial seed point in trachea
    % * assume top of lungs is at end of slice package
    cc = bwconncomp(imclearborder(imclose(Tmask(:,:,end),se)));
    np = cellfun(@length,cc.PixelIdxList);
    slcmask = false(d(1:2));
    slcmask(cc.PixelIdxList{find(np==max(np),1)}) = true;
    slcmask = imdilate(slcmask,se);
    trmask(:,:,end) = slcmask;
    ccs = bwconncomp(slcmask);
    Tmask = Tmask & G;
    
    % Grow from seed mask
    hf = figure;
    hi = imshow(slcmask);
    growing = true;
    slc = d(3)-1;
    while growing
        
        % Find connected regions to previous slice
        seedmask = trmask(:,:,slc+1);
        slcmask = imclose(Tmask(:,:,slc),se);
        cc = bwconncomp(slcmask);
        ind = false(1,cc.NumObjects);
        for i = 1:cc.NumObjects
            if ~any(seedmask(cc.PixelIdxList{i}))
                slcmask(cc.PixelIdxList{i}) = false;
            else
                np = length(cc.PixelIdxList{i});
                nps = 0;
                for j = 1:ccs.NumObjects
                    if any(ismember(cc.PixelIdxList{i},ccs.PixelIdxList{j}))
                        nps = nps + length(ccs.PixelIdxList{j});
                    end
                end
                if np > (A*nps)
                    slcmask(cc.PixelIdxList{i}) = false;
                else
                    ind(i) = true;
                end
            end
        end
        cc.NumObjects = sum(ind);
        cc.PixelIdxList(~ind) = [];
        ccs = cc;
        
        % Set mask slice
        trmask(:,:,slc) = imdilate(slcmask,se);
        set(hi,'CData',trmask(:,:,slc));
        pause(0.01)
        slc = slc - 1;
        
        % Decide whether to stop
        if slc<1 || ~any(slcmask(:))
            growing = false;
        end
    end
    close(hf);
    
    % Remove trachea from lung mask
%     cmiObj.img.mask.morph('erode',[2,2]);
    trmask = cmiObj.img.mask.mat & ~trmask;
    cc = bwconncomp(trmask);
    np = cellfun(@length,cc.PixelIdxList);
    trmask(:) = false;
    trmask(cc.PixelIdxList{np==max(np)}) = true;
    cmiObj.img.mask.merge('replace',trmask);
end