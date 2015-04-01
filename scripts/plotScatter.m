% CMI script
function plotScatter(cmiObj)
%   Input: cmiObj = CMIclass object containing current settings

%% get pixels to put in two columns
img = cmiObj.img.mat(:,:,:,1);
pixelsx = img(cmiObj.img.mask.mat);
img = cmiObj.img.mat(:,:,:,2);
pixelsy = img(cmiObj.img.mask.mat);

pixelsx = pixelsx(1:500:end);
pixelsy = pixelsy(1:500:end);

%% plot colored scatter and labelled axes with emphy and gt cutoffs

%close(1)

ptxmin = -1024;
ptymin = ptxmin;
ptxmax = -500;
ptymax = ptxmax;

h = figure(51)     
dscatter(pixelsx,pixelsy);
hold on
axis([ptxmin ptxmax ptymin ptymax]);
%colormap('bone');
xlabel('Expiration (HU)','FontSize',16);
ylabel('Inspiration (HU)','FontSize',16);
line_gt = [-856 ptymin; -856 ptymax];
line_em = [ptxmin -950; ptxmax -950];
line_ci = [ptxmin (ptxmin-94); ptxmax (ptxmax-94)];

line(line_gt(:,1),line_gt(:,2),'Color','k');
line(line_em(:,1),line_em(:,2),'Color','k');
line(line_ci(:,1),line_ci(:,2),'Color','k');
hold off

