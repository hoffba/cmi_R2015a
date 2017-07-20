% CMIclass function
% Display Max Intensity Projection for current orientation
function showMIP(self,~,~)
hf = figure;
tnc = self.ncolors;
if self.overcheck
    % Background image, set to first image vector
    tclim = self.clim(self.bgvec,:);
    timg = (tnc-1) * ...
        (squeeze(max(self.img.mat(:,:,:,self.bgvec),[],self.orient)) ...
                        - tclim(1)) / diff(tclim);
    timg(timg>tnc) = tnc; % restrict value to color limit max
    image(timg);
    % Overlay Image
    tclim = self.clim(self.vec,:);
    timg = (tnc-1) * ...
        (squeeze(max(self.img.mat(:,:,:,self.vec).*self.img.mask.mat,[],self.orient)) ...
                        - tclim(1)) / diff(tclim);
    timg(timg<1) = 1;
    alpha = 1;
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
image('Parent',gca,'CData',timg+cshift,'AlphaData',adata);
axis image off
set(gca,'YDir','reverse')
tcmap = [eval([self.bgcmap '(' num2str(self.ncolors) ')']);...
    eval([self.cmap '(' num2str(self.ncolors) ')'])];
set(hf,'ColorMap',tcmap,'Name',['MIP: ' self.img.name]);
