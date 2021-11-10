function [prm,info] = CTlung_PRM(img)

%% Load PRM options
info = load(fullfile(fileparts(which('cmi')),'PRMdefs','PRMdef_AllLung_5Color.mat'));

%% Initialize the CMI object
img = ImageClass;
img.setMat(cat(4,exp,insR),{'Exp','Ins'},exp_info.fov,exp_info.orient);
img.mask.merge('replace',mask);
img.prm.setOpts('thresh',prmopt.thresh,...
                'cutoff',prmopt.cutoff,...
                'prmmap',prmopt.prmmap,...
                'cmap',prmopt.cmap,...
                'labels',{'Exp','Ins'},...
                'SPopts',prmopt.SPopts);

%% Generate PRM:
[labels,vals] = img.calcPRM(2);
prm = int8(img.prm.mat);

%% Add results to info:
info.fig = img.prm.hfscatter;
info.regions = labels;
info.pct = vals*100;
