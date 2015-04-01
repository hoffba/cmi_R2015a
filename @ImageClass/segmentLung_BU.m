% ImageClass function
function segmentLung(self,vec)
% Automated lung segmentation algorithm

% Input options
b = self.scaleB(vec);
sr = round(self.dims(1)/100);
dr = 2;
Ti = 550 + b;
Tt = 30 + b;
slc = 174;
prompt = {'Lung Threshold:','Smoothing Radius:',...
    'Trachea Treshold:','Trachea Dilation:','Display:'};
defs = {num2str(Ti),num2str(sr),...
    num2str(Tt),num2str(dr),'on'};
answer = inputdlg(prompt,'Lung Segmentation Options',1,defs);
if ~isempty(answer)
    autoT = strcmpi(answer{1},'auto');
    if ~autoT
        val = str2double(answer{1});
        if ~isempty(val) && ~isnan(val)
            Ti = val;
        end
    end
    val = str2double(answer{2});
    if ~isempty(val) && ~isnan(val)
        sr = round(val);
    end
    val = str2double(answer{3});
    if ~isempty(val) && ~isnan(val)
        Tt = val;
    end
    val = str2double(answer{4});
    if ~isempty(val) && ~isnan(val)
        dr = round(val);
    end
    dopt = strcmpi(answer{5},'on');

    disp('start check'),tic

    tmat = self.mat;
    % Make voxels outside FOV have air HU value
    tmat(tmat<b) = b;

    % First, find optimal threshold
    if autoT
        To = 0; Ti = 50 + b; ct = 0; maxct = 500;
        disp('Finding threshold ...')
        while (Ti~=To) && (ct<=maxct)
            To = Ti;
            ub = mean(tmat(tmat(:)>=To));
            un = mean(tmat(tmat(:)<To));
            Ti = (ub + un)/2;
            ct = ct+1;
            disp(['Count: ' num2str(ct) ' ; Ti: ' num2str(round(Ti))])
        end
    end

    % Find image ends
    tmask = false(self.dims(1:3));
    tmask([1,end],:,:) = true;
    tmask(:,[1,end],:) = true;
    cind = find(tmask);

    % Create initial filtered masks for lung and trachea
    disp('Filtering trachea mask ...')
    th = tic;
    tmask = imfilter(double(tmat<Tt),gauss3D(3,3),0)>0.07;toc(th)
    disp('Filtering lung mask ...')
    th = tic;
    lmask = imerode((tmat<Ti) & ~tmask,ones(3,3,3));toc(th)
    clear tmat

    % Find connected regions
    disp('Finding Regions ...')
    th = tic;
    cc = bwconncomp(lmask);toc(th)
    disp(['Found ' num2str(cc.NumObjects) ' regions ...'])

    % Discard regions that are too small and largest (outside of body)
    n = numel(tmask);
    nvox = cellfun(@numel,cc.PixelIdxList);
    ind = cellfun(@(x)any(ismember(x,cind)),cc.PixelIdxList);
    idx = find(~ind & ((nvox/n)>0.005));
    disp(['Using ' num2str(length(idx)) ' regions ...'])
    lmask2 = false(size(lmask)); lmask3 = lmask2;
    H = fball(7);
    for i = 1:length(idx)
        disp(['Region ' num2str(i) ': ' num2str(nvox(idx(i))) ' voxels (',...
            num2str(nvox(idx(i))*prod(self.voxsz)/10^6) ' L)'])
        th = tic;
        lmask2 = false(size(lmask));
        lmask2(cc.PixelIdxList{idx(i)}) = true;
        disp('dilation ...')
        th2 = tic;
        lmask2 = imdilate(lmask2,H);toc(th2)
        disp('find outside ...')
        th2 = tic;
        cc2 = bwconncomp(~lmask2);toc(th2)
        nvox2 = cellfun(@numel,cc2.PixelIdxList);
        lmask2 = true(size(lmask));
        lmask2(cc2.PixelIdxList{nvox2==max(nvox2)}) = false;
        disp('erosion ...')
        th2 = tic;
        lmask2 = imerode(lmask2,H);toc(th2)
        lmask3 = (lmask3 | lmask2);
        toc(th)
    end
    
    % Find and include inner holes
    %cc = bwconncomp(~lmask3);
    %nvox = cellfun(@numel,cc.PixelIdxList);
    %lmask4 = true(size(lmask));
    %lmask4(cc.PixelIdxList{nvox==max(nvox)}) = false;

    % if dopt
    %     figure;
    %     subplot(),imshow(tmask(:,:,191)
    % end

    self.mask.merge('Replace',lmask3);
    toc
end

