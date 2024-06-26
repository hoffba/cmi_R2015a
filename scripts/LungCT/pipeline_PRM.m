function [prm,info] = pipeline_PRM(exp,exp_info,mask,insR,svname)

%% Load PRM options
info = load(fullfile(fileparts(which('cmi')),'PRMdefs','PRMdef_AllLung_5Color.mat'));

% Check for gapped EXP data
gapchk = nnz(squeeze(std(exp,0,[1,2]))) < size(exp,3);

%% Initialize the CMI object
img = ImageClass;
img.setMat(cat(4,exp,insR),{'Exp','Ins'},exp_info.fov,exp_info.orient,'');
img.mask.merge('replace',mask);
img.prm.setOpts('thresh',info.thresh,...
                'cutoff',info.cutoff,...
                'prmmap',info.prmmap,...
                'cmap',info.cmap,...
                'labels',{'Expiration (HU)','Inspiration (HU)'},...
                'filtchk', true,...
                'filttype', 'median',...
                'filtstr', [3 3 max(3*~gapchk,1)], ...
                'SPopts',info.SPopts,...
                'statchk',false);

%% Generate PRM:
[labels,vals] = img.calcPRM(2);
prm = int8(img.prm.mat);

%% Add title to axes
% tname = exp_info.label;
% if contains(tname,'.')
%     tname = extractBefore(tname,'.');
% end
% title(img.prm.hascatter,tname,'Interpreter','none');

%% Save scatterplot figure
if nargin==5
    pipeline_save_fig(img.prm.hfscatter,svname);
end

%% Add results to info:
info.regions = labels;
info.pct = vals*100;
