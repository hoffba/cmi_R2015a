function [prm,info] = pipeline_PRM(exp,exp_info,mask,insR,svname,gifname)

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
                'filtchk', true,...
                'filttype', 'median',...
                'filtstr', [3 3 3], ...
                'SPopts',info.SPopts);

%% Generate PRM:
[labels,vals] = img.calcPRM(2);
prm = int8(img.prm.mat);

%% Add title to axes
[~,tname] = fileparts(svname);
title(img.prm.hascatter,tname);

%% Save scatterplot figure
print(img.prm.hfscatter,svname,'-dtiff');

%% Add scatterplot to GIF
if nargin==6 && ischar(gifname)
    collate_fig(img.prm.hfscatter,gifname);
end

%% Add results to info:
info.regions = labels;
info.pct = vals*100;
