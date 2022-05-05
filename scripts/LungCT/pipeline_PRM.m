function [prm,info] = pipeline_PRM(exp,exp_info,mask,insR,svname)

%% Load PRM options
info = load(fullfile(fileparts(which('cmi')),'PRMdefs','PRMdef_AllLung_5Color.mat'));

%% Initialize the CMI object
img = ImageClass;
img.setMat(cat(4,exp,insR),{'Exp','Ins'},exp_info.fov,exp_info.orient,'');
img.mask.merge('replace',mask);
img.prm.setOpts('thresh',info.thresh,...
                'cutoff',info.cutoff,...
                'prmmap',info.prmmap,...
                'cmap',info.cmap,...
                'labels',{'Exp','Ins'},...
                'filtchk', 1,...
                'filttype', 'median',...
                'filtstr', [3 3 3], ...
                'SPopts',info.SPopts);

%% Generate PRM:
[labels,vals] = img.calcPRM(2);
prm = int8(img.prm.mat);

%% Save scatterplot figure
print(img.prm.hfscatter,svname,'-dtiff');

%% Add results to info:
info.regions = labels;
info.pct = vals*100;
