function cmiObj0 = fwi_Bruker_proc3(cmiObj0,subjdir)
% Data organized as:
% ~/[subjectID]/[Modality]/[*.mhd]


if nargin==0
    cmiObj0 = CMIclass;
end
if nargin<2
    subjdir = pwd;
end

% Select modality to display:
immod = dir(subjdir);
immod = {immod([immod.isdir]).name};
immod(1:2) = [];
answer = listdlg('ListString',immod);
if isempty(answer)
    return;
else
    immod = immod(answer);
end

% Load display options:
pos = load(fullfile(subjdir,'position.mat'));

% Find background image and VOI:
[~,subjID] = fileparts(subjdir);
bgfname = dir(fullfile(subjdir,'MToff','*.mhd'));
bgfname = fullfile(subjdir,'MToff',bgfname(1).name);
vfname = dir(fullfile(subjdir,'*_VOI.mhd'));
vfname = fullfile(subjdir,vfname(1).name);

% Loop over modality:
tags = {'T2w', 'HighB', 'ADC',  'FatPct', 'R2star', 'MTR' ;
        [],    [],      [0,1],  [2,12],   [0,600],  [0.5,0.8] ; % Color Limits
        [],    [],      0.1,    3,        50,       0.04 };     % PRM threshold
for i = 1:length(immod)
    
    % Load time series modality images:
    j = find(strcmp(immod{i},tags(1,:)),1);
    fname = dir(fullfile(subjdir,immod{i},'*.mhd'));
    fname = fullfile(subjdir,immod{i},{fname.name});
    nf = length(fname);
    bgcheck = ~isempty(tags{3,j});
    if bgcheck
        % Load MToff for background:
        fname = [{bgfname},fname];
    end
    cmiObj0.loadImg(0,fname);
    cmiObj0.loadMask(vfname);
    
    % Set color limits:
    if isempty(tags{2,j})
        clfname = fullfile(subjdir,immod{i},'clim.mat');
        if exist(clfname,'file')
            load(clfname);
        else
            clim = prctile(reshape(cmiObj0.img.mat,prod(cmiObj0.img.dims(1:3)),[]),[25,99]);
            uiwait(msgbox('Adjust contrast.','Adjust Contrast'));
            clim = cmiObj0.clim;
            save(clfname,'clim');
        end
    else
        clim = ones(nf,1)*tags{2,j};
    end
    cmiObj0.setClim([(1:nf)'+bgcheck,clim]);
    
    % Grab VOI mean from unadulterated data:
    vmeans = mean(cmiObj0.img.getMaskVals);
    
    % Adjust image positioning:
    for k = 1:size(pos.rot,1)
        cmiObj0.setView(pos.rot(k,1));
        cmiObj0.imgRotate(pos.rot(k,2));
    end
    cmiObj0.setView(pos.orient);
    cmiObj0.setSlice(pos.slc);
    cmiObj0.imgCrop(1,pos.cropz);
    
    % Interpolate x2
    d = cmiObj0.img.dims(1:3);
    D = d*2;
    fprintf('Interpolating image from (%u %u %u) to (%u %u %u)\n',d,D);
    orient = cmiObj0.orient;
    slc = cmiObj0.slc(orient);
    cmiObj0.imgInterp(D);
    cmiObj0.setView(orient);
    cmiObj0.setSlice(slc*2);
    pause(0.01);
    
    % Adjust VOI:
    cmiObj0.morphMask('close',[3,3]);
    cmiObj0.morphMask('erode',[1,1]);
    cmiObj0.morphMask('open',[1,1]);
    cmiObj0.morphMask('close',[3,3]);
    mask = cmiObj0.img.mask.mat;

end






opts = struct('vdim',{orient},'mdim',{4},'mind',{[]});

% Montage anatomical images (T2w/HighB):
cmiObj0.setOverCheck(false);
cmiObj0.clearMask;
cmiObj0.setCMap('gray');
opts.mind = find(ismember(labi,1:2))';
hf = cmiObj0.genMontage(opts,2);
addLabels(hf,'DWI',[datestrs',{'Date:'}]);
print(hf,fullfile(subjdir,sprintf('%s_DWI_montage.jpg',subjID)),'-djpeg');
cmiObj0.img.mask.merge('replace',mask);

% Generate overlay montages and results:
cmiObj0.setOverCheck(true);
cmiObj0.setCMap('jet');
for i = 3:nt
    opts.mind = find(labi==i)';
    hf = cmiObj0.genMontage(opts,1);
    hf.Children(1).CLim = tags{2,i};
    colormap(hf,'jet');colorbar('FontSize',20);
    pause(0.01);
    addLabels(hf,tags{1,i},...
        [[datestrs';cellfun(@(x)num2str(x,3),num2cell(vmeans(opts.mind)),'UniformOutput',false)],...
        {'Date:';'Mean:'}]);
    print(hf,fullfile(subjdir,sprintf('%s_%s_montageOverlays.jpg',subjID,tags{1,i})),'-djpeg');
end
cmiObj0.setOverCheck(false);
cmiObj0.setCMap('gray');

% Generate PRM montages and scatterplots:
cmiObj0.img.prm.setOpts('cmap',flip(eye(3),2),'filtchk',false,...
    'prmmap',{false(1,2),'PRM_-';[true,false],'PRM_0';true(1,2),'PRM_+'});
thresh = [1,0,1,-1;1,0,1,1];
SPopts = struct('Xvec',1,'Yvec',0,'Xmin',0,'Ymin',1,'Xmax',0,'Ymax',1);
for i = 3:nt
    pth = tags{3,i};
    opts.mind = find(labi==i)';
    thresh(:,1) = opts.mind(1);
    thresh(:,4) = [-1,1]*pth;
    SPopts.Xvec = opts.mind(1);
    SPopts.Xmin = tags{2,i}(1);
    SPopts.Xmax = tags{2,i}(2);
    SPopts.Ymin = tags{2,i}(1);
    SPopts.Ymax = tags{2,i}(2);
    cmiObj0.img.prm.setOpts('thresh',thresh,'SPopts',SPopts);
    
    cmiObj0.activatePRM(true);
    hf = cmiObj0.genMontage(opts,1);
    
    V = cmiObj0.img.getMaskVals(opts.mind);
    [np,ntp] = size(V);
    V = V - repmat(V(:,1),1,ntp);
    V = [ sum(V>=pth,1)/np*100 ; sum(V<-pth,1)/np*100 ];
    addLabels(hf,tags{1,i},...
        [[datestrs';cellfun(@(x)num2str(x,3),num2cell(V),'UniformOutput',false)],...
        {'Date:';'PRM+(%):';'PRM-(%):'}]);
    print(hf,fullfile(subjdir,sprintf('%s_%s_PRMoverlays.jpg',subjID,tags{1,i})),'-djpeg');
end

% Set back to normal modes:
cmiObj0.activatePRM(false);


function addLabels(hf,titleString,labels)
aopts = {'horizontalalignment','center','verticalalignment','middle','fontsize',20,...
    'BackgroundColor','white'};
hf.Name = titleString;
ha = hf.Children(end);
pos1 = ha.Position;
ntp = size(labels,2)-1;
H = 0.08;
Wl = 0.1;
W = pos1(3)/ntp;
posx = pos1(1) + W * ((1:ntp)-1);
for i = 1:ntp
    annotation(hf,'textbox',[posx(i),1-H,W,H],'string',labels(:,i),aopts{:});
end
aopts{2} = 'right';
a = annotation(hf,'textbox',[0,1-H,Wl,H],'string',labels(:,end),'Units','characters',aopts{:});
a.Position(3) = 15;
