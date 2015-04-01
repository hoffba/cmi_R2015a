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

%filt = fspecial('gaussian',[sr,sr],sr);
%se = strel('disk',dr);
se1 = strel('disk',1);
se10 = strel('disk',10);
se9 = strel('disk',9);
% Process each slice separately
d = self.dims(1:3);
np = d(1)*d(2);
tmask = false(d);
h = waitbar(0,'Processing Slices ...');
for i = 1:d(3)
    tm = false(d(1:2));
    tm2 = tm;
    tm([1,end],:) = true;
    tm(:,[1,end]) = true;
    eind = find(tm); % find edge voxels
    tm = tmat(:,:,i)<Ti; % for parenchyma
    tm = imerode(tm,se1);
    cc = bwconncomp(tm);
    nvox = cellfun(@numel,cc.PixelIdxList);
    ind = cellfun(@(x)any(ismember(x,eind)),cc.PixelIdxList);
    idx = find(~ind & (nvox/np)>0.05);
    for j = 1:length(idx)
        ttm = false(size(tm));
        ttm(cc.PixelIdxList{idx(j)}) = true;
        tmask(:,:,i) = (tmask(:,:,i) | imerode(imdilate(ttm,se10),se9));
    end
    waitbar(i/d(3),h)
end
close(h)
toc

self.mask.merge('Replace',tmask);
