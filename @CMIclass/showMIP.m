% CMIclass function
% Display Max Intensity Projection for current orientation
function showMIP(self,~,~)
hf = figure;
tnc = self.ncolors;
oi = 1:3;
oi(oi==self.orient) = [];
voxsz = self.img.voxsz(oi);
if self.overcheck
    % Find indices containing VOI:
    bimg = self.img.mat(:,:,:,self.bgvec);
    oimg = self.img.mat(:,:,:,self.vec).*self.img.mask.mat;
    if self.img.mask.check
        ind = ~squeeze(any(any(self.img.mask.mat,oi(1)),oi(2)));
        switch self.orient
            case 1
                bimg(ind,:,:) = [];
                oimg(ind,:,:) = [];
            case 2
                bimg(:,ind,:) = [];
                oimg(:,ind,:) = [];
            case 3
                bimg(:,:,ind) = [];
                oimg(:,:,ind) = [];
        end
    end
    % Background image, set to first image vector
    tclim = self.clim(self.bgvec,:);
    timg = (tnc-1) * (squeeze(max(bimg,[],self.orient))  - tclim(1)) / diff(tclim);
    timg(timg>tnc) = tnc; % restrict value to color limit max
    image(flipud(timg'));
    % Overlay Image
    tclim = self.clim(self.vec,:);
    timg = (tnc-1) * (squeeze(max(oimg,[],self.orient)) - tclim(1)) / diff(tclim);
    timg(timg<1) = 1;
    alpha = self.dalpha;
    adata = alpha*squeeze(max(self.img.mask.mat,[],self.orient));
    cshift = self.ncolors;
else
    tclim = self.clim(self.vec,:);
    timg = (tnc-1) * (squeeze(max(self.img.mat(:,:,:,self.vec),[],self.orient)) ...
        - tclim(1)) / diff(tclim);
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
