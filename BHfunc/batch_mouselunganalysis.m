
function Cout = batch_mouselunganalysis(C)
% Needs {n x 4} cell array of strings containing mouse lung images/VOIs
%   { Ins , Exp , Exp_VOI , Ins-R , |Jac| }

bpath = '/mnt/cmi/projects/CT_Lung/preclinical/Bleomycin-Fibrosis/Group2/DICOMs';
if nargin==0
    C = cell(20*3,5);
    ct = 1;
    for i = 1:20
        anID = sprintf('Bleo_G2M%02u',i);
        for j = [-1,7,16]
            scanID = sprintf([anID,'_d%02i'],j);
            C(ct,:) = {fullfile(bpath,anID,scanID,[scanID,'_Ins.mhd']),...
                       fullfile(bpath,anID,scanID,[scanID,'_Exp.mhd']),...
                       fullfile(bpath,anID,scanID,[scanID,'_Exp_lungsVOI.mhd']),...
                       fullfile(bpath,anID,scanID,'elxreg','result.mhd'),...
                       fullfile(bpath,anID,scanID,'elxreg','spatialJacobian.mhd')};
            ct = ct+1;
        end
    end
end

% Initialize ImageClass and PRMclass:
img = ImageClass;
img.setPRMopts(1);
pause(0.01);

n = size(C,1);
Cout = [{'bname','Ins-Vol','Ins-Mean','Ins-Med','Exp-Vol','Exp-Mean','Exp-Med',...
         'D','AUC','Ins-Low','Ins-High','Exp-Low','Thresh','Slope'},...
         img.prm.prmmap(:,2)';...
        cell(n,24)];
for i = 1:size(C,1)
    
    [~,bname] = fileparts(fileparts(C{i,3}));
    disp(bname);
    
    if exist(C{i,1},'file') && exist(C{i,2},'file') ...
            && exist(C{i,3},'file') && exist(C{i,4},'file')
        
        tC = cell(1,24);
        
        % Analyze Ins:
%         disp(' - Analyzing Ins ...');
%         img.loadImg(0,C{i,1});
%         mask = quicklungs(img.mat < -200);
%         img.mask.merge('replace',mask&(img.mat<-300));
%         [InsMean,~,InsMed,InsVol,~] = img.getStats(1);
%         tC(2:4) = {InsVol,InsMean,InsMed};

        % Analyze Exp:
        img.loadImg(0,C(i,[2,4,5]));
        img.loadMask(C{i,3});
        img.mask.merge('intersect',img.mat(:,:,:,2)<-200);
        disp(' - Analyzing Ins ...');
        voxvol = prod(img.voxsz);
        tvals = img.getMaskVals([2,3]);
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
%         disp(' - Analyzing Coreg ...');
%         img.loadImg(1,C{i,4});
%         img.loadMask(C{i,3});
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


function omask = quicklungs(mask)
omask = false(size(mask));
omask([1,end],:,:) = true;
omask(:,[1,end],:) = true;
medge = find(omask);
cc = bwconncomp(mask);
i = cellfun(@(x)(~any(ismember(x,medge)))*length(x),cc.PixelIdxList);
i = find(i==max(i),1);
omask(:) = false;
omask(cc.PixelIdxList{i}) = true;

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
