function cmiObj0 = fwi_Bruker_proc(cmiObj0)

% Sort images and set up contrast:
tags = {'T2w', 'HighB', 'ADC',           'FatPct', 'R2star', 'MTR';
        [],    [],      [0,1]*10^-3, [2,12],   [0,400],  [0.5,0.8] };
nt = size(tags,2);
[vmin,vmax] = cmiObj0.img.getColorMinMax;
clim = [vmin,vmax/2];
nvec = cmiObj0.img.dims(4);
labels = cmiObj0.img.labels';
labi = zeros(nvec,1);
for i = 1:nt
    ii = find(~cellfun(@isempty,strfind(labels,strcat('_',tags{1,i}))));
    if ~isempty(ii)
        labi(ii) = i;
        if ~isempty(tags{2,i})
            clim(ii,:) = repmat(tags{2,i},length(ii),1);
        else
            for j = 1:length(ii)
                tvals = reshape(cmiObj0.img.mat(:,:,:,ii(j)),[],1);
                tvals(tvals==0) = [];
                edges = linspace(0,max(tvals)/2,101);
                N = histcounts(tvals,edges);
                NN = gradient(gradient(N));
                mi = max(find(NN==max(NN),1)-5,1);
                tmin = mean(edges([mi,mi+1]));
                clim(ii(j),:) = [tmin,3.5*median(tvals)];
            end
        end
    end
end

[labi,xi] = sort(labi);
clim = clim(xi,:);
cmiObj0.imgOrder(xi);
cmiObj0.setClim([(1:nvec)',clim]);

% Eigen-Transform:
fprintf('Aligning image to axes ...\n');
cmiObj0.eigOrient;
cmiObj0.setView(1);
cmiObj0.imgRotate(-15);
cmiObj0.setView(3);

% Interpolate x2
d = cmiObj0.img.dims(1:3);
D = d*2;
fprintf('Interpolating image from (%u %u %u) to (%u %u %u)\n',d,D);
cmiObj0.imgInterp(D);

% Adjust VOI:
cmiObj0.morphMask('close',[3,3]);
cmiObj0.morphMask('erode',[1,1]);
cmiObj0.morphMask('open',[1,1]);
cmiObj0.morphMask('close',[3,3]);

% Find slice to display:
N = squeeze(sum(sum(cmiObj0.img.mask.mat,1),2));
cmiObj0.setSlice(round(sum( (1:D(3))' .* N / sum(N) )));

% Adjust anatomical images to be similar contrast
uiwait(msgbox('Adjust anatomical contrast.','Contrast'));

% Generate montages and results:
mask = cmiObj0.img.mask.mat;
opts = struct('vdim',{3},'mdim',{4},'mind',{[]});
for i = 1:nt
    opts.mind = find(labi==i)';
    chk = isempty(tags{2,i});
    cmiObj0.setOverCheck(~chk);
    if chk
        cmiObj0.clearMask;
        cmiObj0.setCMap('gray');
    elseif ~cmiObj0.img.mask.check
        cmiObj0.img.mask.merge('replace',mask);
        cmiObj0.setSlice(cmiObj0.slc(cmiObj0.orient));
        cmiObj0.setCMap('jet');
    end
    hf = cmiObj0.genMontage(opts,1);
    hf.Name = tags{1,i};
    if ~chk
        hf.Children(1).CLim = tags{2,i};
        colormap(hf,'jet');colorbar('FontSize',20);
    end
end