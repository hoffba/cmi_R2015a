% CMI script
function mouseLungTest(cmiObj)
%   Input: cmiObj = CMIclass object containing current settings

t = tic;

d = cmiObj.img.dims(1:3);
fov = d.*cmiObj.img.voxsz;
hchk = max(fov)>50;
if hchk % mouse
    Tt = 150;
    Tl = 700;
    r = 5;
    md = 3;
    me = 2;
else % human
    Tt = 50;
    Tl = 700;
    r = 1;
    md = 3;
    me = 2;
end

p = round(0.0002 * prod(d));
img = cmiObj.img.mat(:,:,:,cmiObj.vec);

% Filter image
tmask = false(d); lmask = tmask;
parfor i = 1:d(3)
    if r
        timg = imsharpen(medfilt2(img(:,:,i),[r,r]),'Amount',2,'Radius',r);
    else
        timg = img(:,:,i);
    end
    tmask(:,:,i) = imclearborder(timg < Tt);
    lmask(:,:,i) = imfill(imclearborder(timg < Tl),'holes') & (timg<900);
end
if hchk
    cc = bwconncomp(tmask);
    tmask = false(d);
    np = cellfun(@length,cc.PixelIdxList);
    tmask(cc.PixelIdxList{np==max(np)}) = true;
end
tmask = lmask & ~imdilate(bwareaopen(tmask,p),bwellipsoid(md*ones(1,3)));
tmask = imerode(tmask,bwellipsoid(me*ones(1,3)));
lmask = false(d);
cc = bwconncomp(tmask);
np = cellfun(@length,cc.PixelIdxList);
lmask(cc.PixelIdxList{find(np==max(np),1)}) = true;
cmiObj.img.mask.merge('Replace',lmask);


toc(t)