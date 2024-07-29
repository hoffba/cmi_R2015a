function T = dipg_genstats(num)
try
T = [];

procdir = 'R:\CGalban_Lab\Cancer\DMG_DIPG\20230930_BLau\Processed';
    
if nargin
    id = sprintf('DMG%04d',num);
else
    id = 'DMG';
end
ID = dir(procdir);
ID(~([ID.isdir] & cellfun(@length,{ID.name})==7 & startsWith({ID.name},id) )) = [];
ID = {ID.name}';

% Initialize image object for analysis
img = ImageClass;
img.prm.setOpts('thresh',[2,0,1,-.55; 2,0,1,.55],...
                'prmmap',{logical([0 0]), 'fDM-';...
                          logical([1 0]), 'fDM0';...
                          logical([1 1]), 'fDM+'},...
                'cmap',flip(eye(3)),...
                'statchk',false);

% Loop over subjects
for i = 1:numel(ID)

    regdir = fullfile(procdir,ID{i},'reg');

    % Find ADC maps
    fn_adc = dir(fullfile(regdir,'*.ADC.*nii.gz'));
    fn_adc = {fn_adc.name}';
    caseID0 = extractBefore(fn_adc{1},'.');

    % Find baseline T2w
    fn_T2w = dir(fullfile(regdir,[caseID0,'*.T2w.nii.gz']));

    % Find tumor segmentations
    fn_seg = dir(fullfile(regdir,'*.tumorVOI.nii.gz'));
    fn_seg = {fn_seg.name}';

    % Make sure baseline has required segmentation
    ind = find(startsWith(fn_seg,caseID0),1);
    if ~(isempty(ind) || isempty(fn_T2w))

        % Load baseline T2w for montage background
        img.loadImg(0,fullfile(regdir,fn_T2w(1).name));

        % Load baseline ADC map
        img.loadImg(1,fullfile(regdir,fn_adc{1}));

        % Load baseline mask
        img.loadMask(fullfile(regdir,fn_seg{ind}));
        slcind = sum(img.mask.mat,[1,2]);
        % seg_bl = readNIFTI(fullfile(regdir,fn_seg{ind}));
        % slcind = sum(seg_bl,[1,2]);
        slcind = find(slcind==max(slcind),1);

        % Loop over timepoints
        n = numel(fn_adc);
        M = cell(n+1,1);
        t = table('Size',[n,0]);
        for j = 1:n

            % Find segmentation for this timepoint
            caseID = extractBefore(fn_adc{j},'.');
            ind = find(startsWith(fn_seg,caseID),1);

            if ind
                t.caseID{j} = caseID;
                t.ROI{j} = 'Tumor';

                img.loadImg(1,{fullfile(regdir,fn_adc{j})});
                % img.loadMask(fullfile(regdir,fn_seg{ind}));

                % Grab stats using follow-up tumor volume
                t.Tumor_Vol_cc(j) = prod(img.voxsz)*nnz(img.mask.mat)/10^3;
                vals = img.getMaskVals(3);
                t.Tumor_ADC_mean(j) = mean(vals);
                t.Tumor_ADC_var(j) = var(vals);

                % Stats using intersection of baseline and follow-up tumor segmentation
                % img.mask.merge('intersect',seg_bl);
                % t.fDM_Vol_cc = prod(img.voxsz)*nnz(img.mask.mat)/10^3;
                % vals = img.getMaskVals(3);
                % t.fDM_ADC_mean = mean(vals);
                % t.fDM_ADC_var = var(vals);

                % fDM
                img.prm.setOpts('labels',{['ADC (',caseID0,')'],['ADC (',caseID,')']});
                [labels,vals] = img.calcPRM(3);
                for k = 1:numel(labels)
                    t.(['ADC_',labels{k}])(j) = vals(k)*100;
                end

                % Save PRM scatterplot for QC
                pipeline_save_fig(img.prm.hfscatter,fullfile(regdir,[caseID,'.PRMscatter.tif']));

                % Generate overlay for QC
                ind = any(img.mask.mat,[1,2]);
                dipg_montage(ind,img.mat(:,:,:,1),img.voxsz,img.prm.mat,img.prm.cmap,...
                    fullfile(regdir,[caseID,'.fDM_Overlay.tif']));
                slcIm = dipg_montage(slcind,img.mat(:,:,:,1),img.voxsz,img.prm.mat,img.prm.cmap,...
                    fullfile(regdir,sprintf('%s.Slice%02d.fDM_Overlay.tif',caseID,slcind)));

                % Add date stamp to slice overlay
                d = size(slcIm);
                M{j} = insertText(slcIm,[round(d(2)/2),0],extractAfter(caseID,'_'),...
                    AnchorPoint='CenterTop',FontSize=50,BoxOpacity=0,TextColor='y',Font='Arial Bold');


                % Delete follow-up ADC to prep for next iteration
                img.imgDelete(3);
                % img.mask.clear;
            end
        end

        % Add fDM time plot to montage
        tp = extractAfter(t.caseID,'_');
        tp_dt = datetime(tp,'InputFormat','yyyyMMdd');
        x = days(tp_dt-min(tp_dt));
        figure;
        plot(x,t.("ADC_fDM-"),'-*b',x,100-t.ADC_fDM0,'-*g',x,t.("ADC_fDM+"),'-*r','LineWidth',2);
        ha = gca;
        xlabel(ha,'Days');
        ylabel(ha,'% Tumor Volume');
        legend(ha,'fDM-','100-fDM0','fDM+',"Location",'northwest');
        plotIm = getframe(ha);
        M{n+1} = plotIm.cdata;

        % Montage PRM overlay over time
        imwrite(imtile(M,'GridSize',[round(sqrt((n+1)/2)),NaN]),fullfile(regdir,[caseID,'_TPmontage.tif']));

    end
    % Add stats to table
    T = [T;t];
end

if ~isempty(T)
    writetable(T,fullfile(procdir,'DMG_fDM_Results.xlsx'));
end

catch err
    disp('errchk')
end

function B = dipg_montage(ind,bg,voxsz,lbl,cmap,fname)

lbl = uint8(lbl);
bg_m = imtile(bg(:,:,ind));
lb_m = imtile(lbl(:,:,ind));
B = labeloverlay(bg_m/prctile(bg(:),99.5),uint8(lb_m),"Colormap",cmap,"Transparency",0);
if voxsz(1)~=voxsz(2)
    newd = size(B,[1,2]) .* voxsz(1:2)/min(voxsz(1:2));
    B = imresize(B,newd,'nearest');
end
imwrite(B,fname);





