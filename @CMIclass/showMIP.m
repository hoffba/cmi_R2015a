% CMIclass function
% Display Max Intensity Projection for current orientation
function showMIP(self,x,~)
hf = figure;
tnc = self.ncolors;
oi = 1:3;
oi(oi==self.orient) = [];
voxsz = self.img.voxsz(oi);

fi = 1;
fcnstr = {'max','mean','min'};
fcn = {@(x)squeeze(max(x,[],self.orient)),...
       @(x)squeeze(mean(x,self.orient)),...
       @(x)squeeze(min(x,[],self.orient))};

if nargin>1 && ischar(x) && ismember(x,fcnstr)
    fi = find(strcmp(x,fcnstr));
elseif self.guicheck
    answer = questdlg('What type of projection?','MIP','max','mean','min','max');
    fi = find(strcmp(answer,fcnstr));
end

if self.overcheck
    % Find indices containing VOI:
    bimg = self.img.mat(:,:,:,self.bgvec);
    oimg = self.img.mat(:,:,:,self.vec);
    if self.img.mask.check
        d = self.img.dims(1:3);
        dd = d;
        dd((1:3)~=self.orient) = 1;
        d(self.orient) = 1;
        subs = repmat(reshape(1:prod(d),d),dd);
        subs(~self.img.mask.mat) = prod(d)+1;
        oimg = accumarray(subs(:),oimg(:),[],eval(['@',fcnstr{fi}]));
        oimg = squeeze(reshape(oimg(1:prod(d)),d));
    else
    end
    % Background image, set to first image vector
    tclim = self.clim(self.bgvec,:);
    timg = (tnc-1) * (feval(fcn{fi},bimg)  - tclim(1)) / diff(tclim);
    timg(timg>tnc) = tnc; % restrict value to color limit max
    image(flipud(timg'));
    % Overlay Image
    tclim = self.clim(self.vec,:);
    timg = (tnc-1) * (oimg - tclim(1)) / diff(tclim);
    timg(timg<1) = 1;
    alpha = self.dalpha;
    adata = alpha*squeeze(any(self.img.mask.mat,self.orient));
    cshift = self.ncolors;
else
    timg = feval(fcn{fi},self.img.mat(:,:,:,self.vec));
    tclim = self.clim(self.vec,:);
    timg = (tnc-1) * ( timg - tclim(1) ) / diff(tclim);
    timg(timg>tnc) = tnc;
    adata = 1;
    cshift = 0;
end
image('Parent',gca,'CData',flipud(timg')+cshift,'AlphaData',flipud(adata'));
axis image off
set(gca,'YDir','reverse','DataAspectRatio',[fliplr(voxsz),1])
tcmap = [eval([self.bgcmap '(' num2str(self.ncolors) ')']);...
    eval([self.cmap '(' num2str(self.ncolors) ')'])];
set(hf,'ColorMap',tcmap,'Name',['MIP: ' self.img.name]);

function A = maskedMIP(img,mask,fcn,orient)
fcn = eval(['@',fcn]);
d = size(img);
dd = d;
dd((1:3)~=orient) = 1;
d(orient) = 1;
subs = repmat(reshape(1:prod(d),d),dd);
subs(~mask) = prod(d)+1;
A = accumarray(subs(:),img(:),[],fcn);
A = squeeze(reshape(A(1:prod(d)),d));


