function lmask = segLungMouse(img,opts)

t = tic;

d = size(img);
Tt = 150;
Tl = 700;
r = 5;
md = 3;
me = 2;
if (nargin==2)
    if isfield(opts,'Tt')
        Tt = opts.Tt;
    end
    if isfield(opts,'r')
        r = opts.r;
    end
    if isfield(opts,'Tl')
        Tl = opts.Tl;
    end
    if isfield(opts,'md')
        md = opts.md;
    end
    if isfield(opts,'me')
        me = opts.me;
    end
    clear opts;
end

p = round(0.0002 * prod(d));

% Filter image
tmask = false(d); lmask = tmask;
hw = waitbar(0);
for i = 1:d(3)
    img(:,:,i) = imsharpen(medfilt2(img(:,:,i),[r,r]),'Amount',2,'Radius',r);
    tmask(:,:,i) = imclearborder(img(:,:,i) < Tt);
    lmask(:,:,i) = imfill(imclearborder(img(:,:,i) < Tl),'holes') & (img(:,:,i)<900);
    waitbar(i/d(3),hw)
end
delete(hw)
tmask = lmask & ~imdilate(bwareaopen(tmask,p),bwellipsoid(md*ones(1,3)));
tmask = imerode(tmask,bwellipsoid(me*ones(1,3)));
lmask = false(d);
cc = bwconncomp(tmask);
np = cellfun(@length,cc.PixelIdxList);
lmask(cc.PixelIdxList{find(np==max(np),1)}) = true;

toc(t)



% t = tic;
% 
% % Assumes values are HU + 1000
% Ttrach = 250;
% tdil = 15;
% Tlung = 850;
% lerode = 10;
% dopt = false;
% if (nargin==2)
%     if isfield(opts,'Ttrach')
%         Ttrach = opts.Ttrach;
%     end
%     if isfield(opts,'tdil')
%         tdil = opts.tdil;
%     end
%     if isfield(opts,'Tlung')
%         Tlung = opts.Tlung;
%     end
%     if isfield(opts,'lerode')
%         lerode = opts.lerode;
%     end
%     if isfield(opts,'dopt')
%         dopt = opts.dopt;
%     end
%     clear opts;
% end
% 
% if dopt
%     h = figure;
%     dslc = round(size(img,3)/2);
%     dscale = max(img(:))/2;
% end
% 
% if dopt,figure(h),subplot(2,2,1),imshow(img(:,:,dslc)/dscale),title('Orig');end
% img(img<0) = 0;
% [d(1),d(2),d(3)] = size(img);
% tmask = false(size(img)); % trachea mask
% lmask = tmask; % lungs mask
% 
% % First, smooth image
% p = d(1)/4; % only use middle [p x p] aread of fourier space
% pmask = ones(d(1:2));
% pmask(round(p/2)+1:end-round(p/2),:) = 0;
% pmask(:,round(p/2)+1:end-round(p/2)) = 0;
% for i = 1:d(3)
%     img(:,:,i) = abs(ifft2(fft2(img(:,:,i)).*pmask));
% end
% if dopt,figure(h),subplot(2,2,2),imshow(img(:,:,dslc)/dscale),title('Smooth');end
% 
% % Find trachea:
% disp(' ');
% disp('Finding trachea/bronchi ...');
% tt = tic;
% cc = bwconncomp(img<Ttrach);
% r = inf(1,cc.NumObjects);
% for i = 1:cc.NumObjects
%     [rw,col,slc] = ind2sub(cc.ImageSize,cc.PixelIdxList{i});
%     if any(slc==d(3)) && ~any((rw==1) | (col==1) | (rw==d(1)) | (col==d(2)))
%         r(i) = sqrt((mean(rw)-cc.ImageSize(1)/2)^2 + (mean(col)-cc.ImageSize(2)/2)^2);
%     end
% end
% i = find(r==min(r),1); % closest to center of image
% tmask(cc.PixelIdxList{i}) = true;
% tmask = imdilate(double(tmask),strel('ball',tdil,tdil)) > tdil;
% if dopt,figure(h),subplot(2,2,3),imshow(tmask(:,:,dslc)),title('Bronchi');end
% toc(tt)
% 
% % Find lungs:
% disp('Finding lungs ...')
% tl = tic;
% cc = bwconncomp(img<Tlung);
% r = inf(1,cc.NumObjects);
% for i = 1:cc.NumObjects
%     [rw,col,slc] = ind2sub(cc.ImageSize,cc.PixelIdxList{i});
%     if any(slc==d(3)) && ~any((rw==1) | (col==1) | (rw==d(1)) | (col==d(2)))
%         r(i) = sqrt((mean(rw)-cc.ImageSize(1)/2)^2 + (mean(col)-cc.ImageSize(2)/2)^2);
%     end
% end
% i = find(r==min(r),1); % closest to center of image
% lmask(cc.PixelIdxList{i}) = true;
% lmask = lmask & ~tmask;
% lmask = imerode(double(lmask),strel('ball',lerode,lerode)) > (1-lerode); 
% cc = bwconncomp(lmask);
% np = cellfun(@numel,cc.PixelIdxList);
% lmask = false(size(lmask));
% lmask(cc.PixelIdxList{find(np==max(np),1)}) = true;
% if dopt,figure(h),subplot(2,2,4),imshow(lmask(:,:,dslc)),title('Lungs');end
% toc(tl)
% toc(t)
