
function Cout = batch_mouselunganalysis(C)
% Input:    C = cell array of file names: { Exp-VOI , Exp , Ins-R }
% Output:   Cout = cell array with results

bpath = '/Volumes/projects/CT_Lung_new/preclinical/Bleomycin-Fibrosis/ProcessedData';

% Initialize ImageClass and PRM options:
img = ImageClass;
h = load(fullfile(img.prm.prmdir,'PRMdef_MouseLungALL.mat'));
h.SPopts = struct('Xvec',1,'Yvec',2,...
                  'Xmin',-1000,'Ymin',-1000,...
                  'Xmax',-200,'Ymax',-200,...
                  'show',1,'Nmax',5000);
img.prm.setOpts('thresh',h.thresh,...
                'cutoff',h.cutoff,...
                'cmap',h.cmap,...
                'prmmap',h.prmmap,...
                'SPopts',h.SPopts,...
                'normchk',false);
pause(0.01);

% Initialize output cell array:
n = size(C,1);
% Cout = [{'bname','Ins-Mean','Ins-Med','Exp-Vol','Exp-Mean','Exp-Med',...
%          'D','AUC','Ins-Low','Ins-High','Exp-Low','Thresh','Slope'},...
Cout = [{'bname','Exp-Vol','Ins-Vol','Exp-Mean','Ins-Mean','Exp-Med','Ins-Med',...
         'D','AUC','Exp-Low','Ins-Low','Ins-High','PCA-theta','Thresh','Slope'},...
         img.prm.prmmap(:,2)'];
Cout = [Cout;cell(n,length(Cout))];

for i = 1:n
    
    fpath = fileparts(C{i,1});
    [~,bname] = fileparts(fpath);
    Cout{i+1,1} = bname;
    disp(bname);
    
    if exist(C{i,1},'file') && exist(C{i,2},'file') && exist(C{i,3},'file')
        
        % Find elxreg folder (for SpatialJacobian)
        edir = dir(fullfile(fpath,'elxreg_*'));
        ne = length(edir);
        jchk = false;
        if ne==0
            warning('No elxreg_ folder found for SpatialJacobian.');
        elseif ne>1
            warning('More than one elxreg_ folder found.');
            edir = '';
        else
            edir = fullfile(fpath,edir(1).name);
            jchk = true;
        end
        
        % Load images and mask:
        disp('Loading MHD files ...');
        img.loadImg(0,C(i,2:3));
        img.loadMask(C{i,1});
        if jchk
            img.loadImg(1,fullfile(edir,'spatialJacobian.mhd'));
        end
        
        % Threshold mask to -200HU for analysis:
        img.mask.merge('intersect',max(img.mat(:,:,:,1:2),[],4) < -200);

        % Means:
        voxvol = prod(img.voxsz);
        vals = img.getMaskVals(1:3);
        Cout{i+1,2} = size(vals,1) * voxvol;
        Cout{i+1,4} = mean(vals(:,1));
        Cout{i+1,6} = median(vals(:,1));
        fprintf(' - Exp - Vol= %3.0f ; Mean= %4.0f ; Med= %4.0f\n',Cout{i+1,[2,4,6]});
        if jchk
            Cout{i+1,3} = sum(vals(:,3)) * voxvol;
            Cout{i+1,5} = sum(prod(vals(:,2:3),2))/sum(vals(:,3));
            vals = sortrows(vals,2);
            ind = abs(cumsum(vals(:,3)*voxvol/Cout{i+1,3})-0.5);
            ind = find(ind==min(ind),1);
            Cout{i+1,7} = vals(ind,2);
            fprintf(' - Ins - Vol= %3.0f ; Mean= %4.0f ; Med= %4.0f\n',Cout{i+1,[3,5,7]});
        end
        
        % PP-plot
        clearText(findobj('Name','Histogram Fit (X)'));
        clearText(findobj('Name','Histogram Fit (Y)'));
        s = img.ppPlot(2,[ 0 1 1 1 1 ]);
        Cout(i+1,8:13) = {s.D,s.AUC,s.xEsts(4),s.yEsts(4),s.yEsts(5),s.PCA.theta};
        pause(0.01);
        
        % CI
        s = img.calcPRMCI(2);
        Cout(i+1,14:15) = {abs(s.CI95),s.slope};
        pause(0.01);
        
        % PRM
        [~,PRMvals] = img.calcPRM(2);
        Cout(i+1,16:end) = num2cell(PRMvals)';
        pause(0.01);

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
