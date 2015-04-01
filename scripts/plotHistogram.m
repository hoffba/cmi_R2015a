% CMI script
function plotHistogram(cmiObj)
%   Input: cmiObj = CMIclass object containing current settings

%% get pixels to put in two columns
img = cmiObj.img.mat(:,:,:,1);
pixelsx = img(cmiObj.img.mask.mat);
img = cmiObj.img.mat(:,:,:,2);
pixelsy = img(cmiObj.img.mask.mat);


%% plot colored scatter and labelled axes with emphy and gt cutoffs

%close(1)

ptxmin = -1024;
ptymin = 0;
ptxmax = -500;
ptymax = 2;

h = figure(52)     
[xnelements xcenters] = hist(pixelsx,150);
plot(xcenters,xnelements,'Color','m');
hold on
[ynelements ycenters] = hist(pixelsy,150);
plot(ycenters,ynelements,'Color','c');

% Set additional plot characteristics
xlim([ptxmin ptxmax]);
xlabel('Image (HU)','FontSize',16);
ylabel('Number of Voxels','FontSize',16);
set(gca,'XTick',-1000:100:-500);
%set(gca,'YTick',0:1:5); % This is arbitrary and only for right now!

hold off
