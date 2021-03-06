function cmiObj0 = fwi_Bruker_proc2(cmiObj0,subjdir)

if nargin==0
    cmiObj0 = CMIclass;
end
if nargin<2
    subjdir = pwd;
end

% Load images:
[~,subjID] = fileparts(subjdir);
[stat,datestrs] = loadMFdata(cmiObj0,subjdir);

% Sort images and set up contrast:
tags = {'T2w', 'HighB', 'ADC',  'FatPct', 'R2star', 'MTR' ;
        [],    [],      [0,1],  [2,12],   [0,600],  [0.5,0.8] ;
        [],    [],      0.1,   3,        50,       0.04 };
nt = size(tags,2);
[vmin,vmax] = cmiObj0.img.getColorMinMax;
clim = [vmin,vmax/2];
nvec = cmiObj0.img.dims(4);
labels = cmiObj0.img.labels';
labi = zeros(nvec,1);
for i = 1:nt
    ii = find(~cellfun(@isempty,strfind(labels,tags{1,i})));
    if ~isempty(ii)
        labi(ii) = i;
        if ~isempty(tags{2,i})
            clim(ii,:) = repmat(tags{2,i},length(ii),1);
        end
    end
end

[labi,xi] = sort(labi);
clim = clim(xi,:);
cmiObj0.imgOrder(xi);
cmiObj0.setClim([(1:nvec)',clim]);
vmeans = mean(cmiObj0.img.getMaskVals);

% Set color limits for DWI:
ntp = length(datestrs);
fname = fullfile(subjdir,'clim.mat');
if exist(fname,'file')
    t = load(fname,'clim');
    tclim = t.clim;
else
    tclim = zeros(0,2);
end
ntpt = round(size(tclim,1)/2);
itp = (ntpt+1:ntp)+1;
clim((1:ntpt)+1,:) = tclim(1:ntpt,:);
clim((1:ntpt)+1+ntp,:) = tclim((1:ntpt)+ntpt,:);
cmiObj0.setClim([(1:nvec)',clim]);
cmiObj0.h.slider_4D.Enable = 'off'; pause(0.01);
for i = 1:length(itp)
    cmiObj0.setVec(itp(i));
    uiwait(msgbox('Adjust contrast.','Adjust Contrast'));
end
cmiObj0.h.slider_4D.Enable = 'on'; pause(0.01);
clim = cmiObj0.clim(1+(1:ntp*2),:);
save(fullfile(subjdir,'clim.mat'),'clim');

% Load position parameters previously saved by user:
cmiObj0.setVec(1);
if exist(fullfile(subjdir,'position.mat'),'file')
    pos = load(fullfile(subjdir,'position.mat'));
    if ~all(ismember({'rot','orient','slc','cropz'},fieldnames(pos)))
        error('Missing position information.');
    end
    for i = 1:size(pos.rot,1)
        cmiObj0.setView(pos.rot(i,1));
        cmiObj0.imgRotate(pos.rot(i,2));
    end
    cmiObj0.setView(pos.orient);
    cmiObj0.setSlice(pos.slc);
    ind = 1:3; ind(pos.orient) = [];
    cmiObj0.imgCrop(ind(2),pos.cropz);
else
    error('Position info not saved.');
end

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

% mask = cmiObj0.img.mask.mat;
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
if ntp>1
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
end



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
