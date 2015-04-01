function lmask = growMouseLungsVOI(img,slc,slcmask,T,Tl)
% Removes airways from existing lung mask
% Assumes a filtered image is loaded (low noise from NOVA filter)

lmask = [];
d = size(img);
if length(d)<4
    % Set growing parameters:
    if nargin<2
        T = -740; % airways threshold
    end
    if nargin<3
        Tl = -200; % lung threshold
    end
    A = 1.5; % 2D region size threshold
    se = strel('disk',5);
    Tmask = img(:,:,:,1) < T; % Exp image
    
    % Initialize 3D matrices:
    trmask = false(d);
    
    % Find initial seed point in trachea
    % * assume top of lungs is at beginning of slice package
%     cc = bwconncomp(imclearborder(imclose(Tmask(:,:,slc-1),se)));
%     stats = regionprops(cc,{'Centroid'});
%     C = d(1:2)/2;
%     R = cellfun(@(x) sqrt(sum((x-C).^2)),{stats(:).Centroid});
%     slcmask = false(d(1:2));
%     slcmask(cc.PixelIdxList{find(R==min(R),1)}) = true;
    slcmask = imdilate(slcmask,se);
    trmask(:,:,slc-1) = slcmask;
    ccs = bwconncomp(slcmask);
    
    % Grow from seed mask
    growing = true;
    while growing
        
        % Find connected regions to previous slice
        seedmask = trmask(:,:,slc-1);
        slcmask = Tmask(:,:,slc);
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
        slc = slc + 1;
        
        % Decide whether to stop
        if slc<1 || ~any(slcmask(:))
            growing = false;
        end
    end
    
    se = bwellipsoid([3,3,3]);
    lmask = imerode(img<Tl,se);
    lmask(:,:,end) = zeros(d(1:2));
    lmask(:,:,1) = zeros(d(1:2));
    lmask = imclearborder(lmask) & ~trmask;

    ii = squeeze(any(any(lmask,1),2));
    imin = find(ii,1); imax = find(ii,1,'last');
    parfor i = imin:imax
        lmask(:,:,i) = medfilt2(lmask(:,:,i),[10,10]);
    end
    cc = bwconncomp(lmask);
    np = cellfun(@length,cc.PixelIdxList);
    lmask(:) = false;
    lmask(cc.PixelIdxList{np==max(np)}) = true;
end