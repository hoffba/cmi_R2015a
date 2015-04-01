% CMI script
function plotScatterGYR(cmiObj)
%   Input: cmiObj = CMIclass object containing current settings
% 
%   Description: Assumes the cmi program has at least two images and a VOI
%   loaded. Creates a scatter plot with CT lung PRM colored background.
%   Sets limits on plot region to lung values.
%
%   Revision history:
%   16 December 2013 jlb created 
%   26 Feb 2014 jlb added dialogue for how much to subsample in order to
%                   get smaller eps plots

%% check for existence of 2 images and mask
% if ~cmiObj.img.mat.check
%     disp('Error: need a mask');
%     return;
% elseif cmiObj.img.dims(1,3)<2
%     disp('Error: need at least two images');
%     return;
% end

%% SET EVERY Xth Pixel (default is 500)
pixel_select_every = 500;
prompt = {'Select every N voxels for plot:'};
dlg_title = 'Input';
num_lines = 1;
def = {'500'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
if ~isempty(answer)
    pixel_select_every = str2num(answer{:});
end

%% get pixels to put in two columns
img = cmiObj.img.mat(:,:,:,1);
pixelsx = img(cmiObj.img.mask.mat);
img = cmiObj.img.mat(:,:,:,2);
pixelsy = img(cmiObj.img.mask.mat);

pixelsx = pixelsx(1:pixel_select_every:end);
pixelsy = pixelsy(1:pixel_select_every:end);

%% plot colored scatter and labelled axes with emphy and gt cutoffs

%close(1)

ptxmin = -1024;
ptymin = ptxmin;
ptxmax = -500;
ptymax = ptxmax;
xsize = ptxmax-ptxmin+1;
ysize = ptymax-ptymin+1;
ythresh = 1024-950+1;
xthresh = 1024-856+1;

h = figure(51) 

% create colored background for plot
img = zeros(ptxmax-ptxmin+1,ptymax-ptymin+1);
cimg = zeros(size(img,1),size(img,2),3);
cimg(:,:,:) = 1;
cimg(1:ythresh,1:xthresh,:) = 0.2; %R
cimg(1:ythresh,1:xthresh,1) = 1; %R
cimg(ythresh:ysize,1:xthresh,:) = 1; %Y
cimg(ythresh:ysize,1:xthresh,3) = 0.2; %Y
cimg(ythresh:ysize,xthresh:xsize,:) = 0.2; %G
cimg(ythresh:ysize,xthresh:xsize,2) = 1; %G
image(cimg,'XData',-1024,'YData',-1024); % This works, but covers up graph!!!!
set(gca,'ydir','normal'); % for image display fixing
hold on

% scatter plot
dscatter(pixelsx,pixelsy);
hold on

% axis limits and labels
% Note: 48 & 32 is too big to downsize to about 1"
% Note: 18 and 16 is too small ""
% Note: 32 and 18 not too bad (18 a little small, so is 32)
axis([ptxmin ptxmax ptymin ptymax]);
xlabel('Expiration (HU)','FontSize',36);
ylabel('Inspiration (HU)','FontSize',36);
set(gca,'XTick',[-1000, -856, -500],'FontSize',24);
set(gca,'YTick',[-1000, -950, -500],'FontSize',24);

% cutoff lines
line_gt = [-856 ptymin; -856 ptymax];
line_em = [ptxmin -950; ptxmax -950];
line_ci = [ptxmin (ptxmin-94); ptxmax (ptxmax-94)];
line(line_gt(:,1),line_gt(:,2),'Color','k');
line(line_em(:,1),line_em(:,2),'Color','k');
%line(line_ci(:,1),line_ci(:,2),'Color','k');



%image(img);

hold off;




