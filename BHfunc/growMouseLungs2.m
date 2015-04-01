function growMouseLungs(cmiObj,T)
% Removes airways from existing lung mask
% Assumes a filtered image is loaded (low noise from NOVA filter)

if cmiObj.img.check
    t = tic;
    % Set growing parameters:
    if nargin<2
        T = -740; % Image threshold
    end
    A = 1.3; % 2D region size threshold
    se = strel('disk',5);
    Tmask = cmiObj.img.mat(:,:,:,1) < T; % Exp image
    
    % Initialize 3D matrices:
    d = cmiObj.img.dims(1:3);
    slc = d(3)-11;
    trmask = false(d);
    
    % Find initial seed point in trachea
    % * assume top of lungs is at end of slice package
    cc = bwconncomp(imclearborder(imclose(Tmask(:,:,slc+1),se)));
    stats = regionprops(cc,{'Centroid'});
    C = d(1:2)/2;
    R = cellfun(@(x) sqrt(sum((x-C).^2)),{stats(:).Centroid});
%     np = cellfun(@length,cc.PixelIdxList);
    slcmask = false(d(1:2));
%     slcmask(cc.PixelIdxList{find(np==max(np),1)}) = true;
    slcmask(cc.PixelIdxList{find(R==min(R),1)}) = true;
    slcmask = imdilate(slcmask,se);
    trmask(:,:,slc+1) = slcmask;
    ccs = bwconncomp(slcmask);
    
    % Grow from seed mask
    hf = figure;
    hi = imshow(slcmask);
    growing = true;
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
    cmiObj.img.mask.merge('replace',trmask);
    answer = questdlg('Continue (Check trachea VOI)?');
    
    if strcmp(answer,'Yes')
        % Remove trachea from lung mask
        hf = waitbar(0,'Thresholding for lungs ...');
        se = bwellipsoid([3,3,3]);
        lmask = imerode(cmiObj.img.mat(:,:,:,1)<-200,se);
        lmask(:,:,end) = zeros(d(1:2));
        lmask = imclearborder(lmask) & ~trmask;
        waitbar(0.5,hf,'Extracting lung regions ...');
        cc = bwconncomp(lmask);
        np = cellfun(@length,cc.PixelIdxList);
        lmask(:) = false;
        lmask(cc.PixelIdxList{np==max(np)}) = true;
        waitbar(1,hf,'Filtering lung mask ...');
        ii = squeeze(any(any(lmask,1),2));
        imin = find(ii,1); imax = find(ii,1,'last');
        if (exist('matlabpool')==2) && (matlabpool('size')==0)
            matlabpool;
        end
        parfor i = imin:imax
            lmask(:,:,i) = medfilt2(lmask(:,:,i),[10,10]);
        end
        delete(hf);
        cmiObj.img.mask.merge('replace',lmask);
        toc(t)
    end
end