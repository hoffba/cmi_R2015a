
function C = batch_LungPRM_Mosaic
% Processes lung image data for the Mosaic clinical project
% Requires: Exp/Ins, VOIs for both, Coregistered Ins, SpatialJacobian
%   Coregister using "regLungs.m"

% Hard-coded base path to finde processed MHD files:
bpath = '/mnt/cmi/projects/CT_Lung_new/clinical/data/UM/Mosaic Project-17Nov2014/processed/';
dnames = dir(bpath);
dnames = {dnames([dnames(:).isdir]).name};

% Initialize ImageClass and PRMclass:
img = ImageClass;
load(fullfile(fileparts(which('cmi')),'PRMdefs','PRMdef_AllLung.mat'))
img.prm.setOpts('thresh',thresh,...
                'cutoff',cutoff,...
                'cmap',cmap,...
                'prmmap',prmmap,...
                'labels',{'Exp','Ins'},...
                'npmaxscat',5000,...
                'prmscatter',true,...
                'normchk',true);
pause(0.01);

% Initialize output cell array:
Head = [{'bname','Ins-Vol','Ins-Mean','Ins-Med','Exp-Vol','Exp-Mean','Exp-Med',...
         'D','AUC','Ins-Low','Ins-High','Exp-Low','Thresh','Slope'},...
         img.prm.prmmap(:,2)'];
     
% Loop over all patient directories:
C = cell(0,length(Head));
for i = 1:length(dnames)
    
    bname = dnames{i};
    disp(bname);
    
    % Find images:
    go = true;
    inames = dir(fullfile(bpath,bname,'*.mhd'));
    inames = {inames(:).name};
    oexpfn = find(~cellfun(@isempty,strfind(inames,'_Exp.mhd')),1);
    if isempty(oexpfn)
        warning('*_Exp.mhd not found.');
        go = false;
    else
        oexpfn = fullfile(bpath,bname,inames{oexpfn});
    end
    oinsfn = find(~cellfun(@isempty,strfind(inames,'_Ins.mhd')),1);
    if isempty(oinsfn)
        warning('*_Ins.mhd not found.');
        go = false;
    else
        oinsfn = fullfile(bpath,bname,inames{oinsfn});
    end
    expvoifn = find(~cellfun(@isempty,strfind(inames,'_Exp_VOI.mhd')),1);
    if isempty(expvoifn)
        warning('*_Exp_VOI.mhd not found.');
        go = false;
    else
        expvoifn = fullfile(bpath,bname,inames{expvoifn});
    end
    insvoifn = find(~cellfun(@isempty,strfind(inames,'_Ins_VOI.mhd')),1);
    if isempty(insvoifn)
        warning('*_Ins_VOI.mhd not found.');
        go = false;
    else
        insvoifn = fullfile(bpath,bname,inames{insvoifn});
    end
    [~,elxdir] = fileparts(oinsfn);
    elxdir = ['elxreg_',elxdir];
    rinsfn = fullfile(bpath,bname,elxdir,'result.mhd');
    if ~exist(rinsfn,'file')
        warning('Coregistered Ins not found.');
        go = false;
    end
    jacfn = fullfile(bpath,bname,elxdir,'spatialJacobian.mhd');
    if ~exist(jacfn,'file')
        warning('Spatial Jacobian not found.');
        go = false;
    end
    
    if go 
        
        % Analyze Ins:
        disp(' - Analyzing Ins ...');
        img.loadImg(0,oinsfn);
        img.loadMask(insvoifn);
        [InsMean,~,InsMed,InsVol,~] = img.getStats(1);
        fprintf(' - - Vol= %3.0f ; Mean= %4.0f ; Med= %4.0f\n',InsVol,InsMean,InsMed);
        C(i,2:4) = {InsVol,InsMean,InsMed};

        % Analyze Exp:
        disp(' - Analyzing Exp ...');
        img.loadImg(0,{oexpfn,rinsfn});
        img.loadMask(expvoifn);
        [ExpMean,~,ExpMed,ExpVol,~] = img.getStats(1);
        fprintf(' - - Vol= %3.0f ; Mean= %4.0f ; Med= %4.0f\n',ExpVol,ExpMean,ExpMed);
        C(i,5:7) = {ExpVol,ExpMean,ExpMed};

        % Analyze Coreg:
        disp(' - Analyzing Coreg ...');
        clearText(findobj('Name','Histogram Fit (X)'));
        clearText(findobj('Name','Histogram Fit (Y)'));
        PPstats = img.ppPlot(2);
        C(i,8:12) = {PPstats.D,PPstats.AUC,PPstats.yEsts(4),PPstats.yEsts(5),PPstats.xEsts(4)};
        pause(0.01);
        CIstats = img.calcPRMCI(2);
        C(i,13:14) = {abs(CIstats.CI95),CIstats.slope};
        pause(0.01);
        [~,PRMvals] = img.calcPRM(2);
        C(i,15:end) = num2cell(PRMvals)';
        pause(0.01);

    else
        disp('Skipping patient due to missing data.');
    end
end
C = [Head;C];
cmi_csvwrite(fullfile(bpath,'batch.csv'),C);

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
