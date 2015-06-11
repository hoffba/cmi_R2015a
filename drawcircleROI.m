%% Initialize series
count = 0;
ellipseRadPixel = 10;

%% Initialize new case

% You should: 
% (1) Read in case
% (2) Scale down HU if needed
% (3) THEN execute this section to initialize for new case

% Select cmiObj image figure and get information
figure(cmiObj0.hfig);
h = gca;
pos = get(h,'position');
xlim = get(h,'xlim');
ylim = get(h,'ylim');
dim = cmiObj0.img.dims(1:2);

% Center for initial roi placement
cntr = [(xlim(2)-xlim(1))/2 (ylim(2)-ylim(1))/2];
cntr_in_pixels = [dim(1)*((cntr(1)-xlim(1))/(xlim(2)-xlim(1)))...
   dim(2)*((cntr(2)-ylim(1))/(ylim(2)-ylim(1)))] ;
pos = [cntr ellipseRadPixel ellipseRadPixel]; % initial upper left corner & bounding rectangle
count = count+1;



%% Cycle thru specific case and find ROI's

% initialiize
roiType{1} = 'outside_air';
roiType{2} = 'trachea_air';
roiType{3} = 'aorta_blood';

% cycle thru more than one study
roiName{count} = cmiObj0.img.name;

for i=1:3
    
    % Create roi
    hellipse = imellipse(h,pos);
    if i==3; hellipse.setColor('r'); elseif i==1; hellipse.setColor([0.5 0.5 1]); end;
    fprintf('Move ROI to cover %s (double click when placed)\n',roiType{i});
    roiplace = wait(hellipse);
    fprintf('DONE: %s placement\n',roiType{i});
    roipos = hellipse.getPosition();
    
    % Get roi stats
    
    % convert axes coordinates to image pixel coordinates
    % I am also rounding to integers (ncol x, nrow y, XData, axes coord)
    %roiloc = hellipse.getPosition();
    %p = roiloc(1:2);
    pMat = [dim(1)*((roipos(1)-xlim(1))/(xlim(2)-xlim(1)))...
        dim(2)*((roipos(2)-ylim(1))/(ylim(2)-ylim(1)))] ;
    pscale = [(xlim(2)-xlim(1))/dim(1) (ylim(2)-ylim(1))/dim(2)];
    
    pdMat = pscale.*roipos(3:4);
    
    % get mask and stats
    pcMat = pMat + pdMat;
    [rr cc] = meshgrid(1:dim(1),1:dim(2));
    C = sqrt((rr-pcMat(1)).^2+(cc-pcMat(2)).^2)<=pdMat(1);
%     figure(4);
%     imshow(C);
    currImg = cmiObj0.img.mat(:,:,cmiObj0.slc(3));
    mHU = mean(currImg(C));
    fprintf('DONE with Roi %s mean HU=%f\n',roiType{i},mHU);
    
    
    roiData(count,i) = mHU;
    
    % delete roi
    hellipse.delete();

end

%% Write out data

% Output file for .csv prm info
[out_tname, out_tpath] = uiputfile('*.csv','Select output .csv file');
if isequal(out_tpath,0) || isequal(out_tname,0)
    fprintf('Error: no .csv file name selected. Exiting script.\n');
    return;
end
% Open file
csvfname= fullfile(out_tpath,out_tname);
fid = fopen(csvfname,'w');
if fid==-1
    fprintf('****ERROR opening %s\n',out_tname);
    return;
end

% output ['Patient' outlabels] and [patient_id outvals]
fprintf(fid,'Patient,%s,%s,%s\n',roiType{1},roiType{2},roiType{3});
for tt=1:size(roiData,1)
    fprintf(fid,'%s,%f,%f,%f\n',roiName{tt},roiData(tt,1),roiData(tt,2),roiData(tt,3));
end % filenames
fclose(fid);


%% IGNORE
% this doesn't work because the fig is not image data
% [nrows,ncols] = size(get(h,'CData'));
% xdata = get(h,'XData');
% ydata = get(h,'YData');
% px = axes2pix(ncols,xdata,30);
% py = axes2pix(nrows,ydata,30);   
        
% this doesn't work because h is not an image type
%mellipse = hellipse.createMask(h);