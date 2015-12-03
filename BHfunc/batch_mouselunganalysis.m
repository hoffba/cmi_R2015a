
function Cout = batch_mouselunganalysis(C)
% Input:    C = cell array of file names: { Exp-VOI , Exp , Ins-R }
% Output:   Cout = cell array with results


% Initialize ImageClass and PRM options:
img = ImageClass;
h = load(fullfile(img.prm.prmdir,'PRMdef_MouseLungALL.mat'));
img.prm.setOpts('thresh',h.thresh,...
                'cutoff',h.cutoff,...
                'cmap',h.cmap,...
                'prmmap',h.prmmap,...
                'SPopts',h.SPopts);
pause(0.01);

% Initialize output cell array:
n = size(C,1);
Cout = [{'bname','Ins-Mean','Ins-Med','Exp-Vol','Exp-Mean','Exp-Med',...
         'D','AUC','Ins-Low','Ins-High','Exp-Low','Thresh','Slope'},...
         img.prm.prmmap(:,2)'];
Cout = [Cout;cell(n,length(Cout))];

for i = 1:n
    
    bname = strsplit(C{i,1},filesep);
    bname = bname{end-1};
    disp(bname);
    
    if exist(C{i,1},'file') && exist(C{i,2},'file') && exist(C{i,3},'file')
        
        tC = cell(1,24);
        
        % Load images and mask:
        img.loadImg(0,C(i,2:3));
        img.loadMask(C{i,1});
        
        % Threshold mask to -200HU for analysis:
        img.mask.merge('intersect',prod(img.mat<-200,4));

        % Analyze Exp:
        disp(' - Analyzing Ins ...');
        voxvol = prod(img.voxsz);
        tvals = img.getMaskVals(1:2);
        InsVol = sum(tvals(:,2)) * voxvol;
        InsMean = sum(prod(tvals,2)) / InsVol * voxvol;
        % Compute weighted median:
        tvals = sortrows(tvals,1);
        ind = abs(cumsum(tvals(:,2)/InsVol*voxvol)-0.5);
        ind = find(ind==min(ind),1);
        InsMed = tvals(ind,1);
        fprintf(' - - Vol= %3.0f ; Mean= %4.0f ; Med= %4.0f\n',InsVol,InsMean,InsMed);
        tC(2:4) = {InsVol,InsMean,InsMed};
        disp(' - Analyzing Exp ...');
        [ExpMean,~,ExpMed,ExpVol,~] = img.getStats(1);
        fprintf(' - - Vol= %3.0f ; Mean= %4.0f ; Med= %4.0f\n',ExpVol,ExpMean,ExpMed);
        tC(5:7) = {ExpVol,ExpMean,ExpMed};

        % Analyze Coreg:
        clearText(findobj('Name','Histogram Fit (X)'));
        clearText(findobj('Name','Histogram Fit (Y)'));
        PPstats = img.ppPlot(2);
        tC(8:12) = {PPstats.D,PPstats.AUC,PPstats.yEsts(4),PPstats.yEsts(5),PPstats.xEsts(4)};
        pause(0.01);
        CIstats = img.calcPRMCI(2);
        tC(13:14) = {abs(CIstats.CI95),CIstats.slope};
        pause(0.01);
        [~,PRMvals] = img.calcPRM(2);
        tC(15:end) = num2cell(PRMvals)';
        pause(0.01);

        % Add stats to Cout:
        Cout(i+1,:) = tC;
    else
        disp(' - Doesn''t exist');
    end
end
cmi_csvwrite(fullfile(bpath,'batch.csv'),Cout);

function clearText(h)
if ~isempty(h) && ishandle(h)
    h = allchild(h);
    h = h(strcmp('axes',get(h,'Type')));
    for i = 1:length(h)
        th = allchild(h(i));
        if all(strcmp('text',get(th,'Type')))
            delete(th);
        end
    end
end
