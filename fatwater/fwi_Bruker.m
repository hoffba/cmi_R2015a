function fwi_Bruker(studydir)
% Assumed acquisition:



interpchk = true;
interpm = 'cubic';

[~,studystr] = fileparts(studydir);
studystr = strsplit(studystr,'_');
fnout = strjoin(studystr([3:end-2,1]),'_');
svdir = fullfile(studydir,'MHD');
if ~isdir(svdir)
    mkdir(svdir);
end

% Catalog study directory:
[C,hfig] = mriCatalog(studydir,1);

% Find relevant images:
ind = {find(cellfun(@(x)~isempty(strfind(x,'_Fat')),C(:,2)))',...
        find(cellfun(@(x)~isempty(strfind(x,'_DtiStandard')),C(:,2)))',...
        find(cellfun(@(x)~isempty(strfind(x,'_MT')),C(:,2)))'};
answer = inputdlg({'Fat:','Diffusion:','MT'},'MF processing',1,cellfun(@num2str,ind,'UniformOutput',false));
delete(hfig);

%% Fat:
% 3D coronal
ser = C(str2double(strsplit(answer{1})),1);
for i = 1:length(ser)
    fprintf('Loading series: %s\n',ser{i});
    [timg,~,tfov,p] = readBrukerMRI(fullfile(studydir,ser{i},'fid'));
    [td(1),td(2),td(3),~] = size(timg);
    if i==1
        I = coilCombine(timg);
        fov = tfov([2,1,3]);
        te = p.EffectiveTE;
        f0 = p.PVM_FrqRef(1);
        d = td;
    elseif all(d==td)
        I = cat(5,I,coilCombine(timg));
        te = cat(2,te,p.EffectiveTE);
    else
        fprintf('Dimensions mismatch - Series %s - ( %u %u %u ) ~= Orig( %u %u %u )\n',ser{i},td,d);
    end
end

% Sort images:
[te,xi] = sort(te);
I = I(:,:,:,:,xi);

% Calculate Fat:
% Fat model:
gyro = 42.57747892; % MHz/T
B0 = f0/gyro;
aF = [0.087 0.693 0.128 0.004 0.039 0.048]; % Fat amplitudes
wF = [-3.8, -3.4, -2.67, -1.93, -0.4, 0.60]; % rad/ms
% Diego Hernando
% Perform 2D complex analysis using Graph Cuts
fprintf('Computing fat using 2D complex graph cuts ...\n');
tdata = struct('images',[],... % [nx,ny,nz,nC,nE,nAcq]
               'TE',te/1000,... % need to be in seconds
               'FieldStrength',B0,...
               'PrecessionIsClockwise',0);
params = struct('species',struct('name',{'water','fat'},'frequency',{0,wF},'relAmps',{1,aF}),...
                'size_clique',1,...
                'range_r2star',[0 2500],...
                'NUM_R2STARS',101,...
                'range_fm',[-1000 1000],...
                'NUM_FMS',101,...
                'NUM_ITERS',40,...
                'SUBSAMPLE',4,...
                'DO_OT',1,...
                'LMAP_POWER',2,...
                'lambda',0.05,...
                'LMAP_EXTRA',0.05,...
                'TRY_PERIODIC_RESIDUAL',0);
% 2D slice by slice:
di = 3;
dord = [setdiff(1:3,di),di];
I = permute(I,[dord,4,5]);
F = zeros(d(dord)); W = F; R2s = F; fmap = F;
for i = 1:d(di)
    fprintf('  Slice %u / %u\n',i,d(di));
    tdata.images = I(:,:,i,:,:);
    results = fw_i2cm1i_3pluspoint_hernando_graphcut(tdata,params);
    % Collect results: (F/W are swapped for some reason)
    F(:,:,i) = results.species(1).amps;
    W(:,:,i) = results.species(2).amps;
    R2s(:,:,i) = results.r2starmap;
    fmap(:,:,i) = results.fieldmap;
end
[~,dord] = sort(dord);
I = permute(I,[dord,4,5]);
W = abs(permute(W,dord));
F = abs(permute(F,dord));
R2s = permute(R2s,dord);
fmap = permute(fmap,dord);
FF = computeFF(W,F);
% Save results:
save(fullfile(svdir,sprintf('%s_Fat.mat',fnout)),'I','te','fov','params');
saveMHD(fullfile(svdir,sprintf('%s.mhd',fnout)),...
    flip(cat(4,FF,abs(W),abs(F),R2s,fmap,mean(abs(I),5)),1),...
    strcat('FWI_',{'FatPct','Water','Fat','R2star','FieldMap','MGEmean'}),fov);

%% Determine interpolation image geometry:
% Assumes all images use the same offsets
if interpchk
    fov0 = fov;
    d0 = d;
    v = fov./d;
    vd = (d-1)/2 .* v;
    [Xq,Yq,Zq] = meshgrid(-vd(2):v(2):vd(2),-vd(1):v(1):vd(1),-vd(3):v(3):vd(3));
end

%% Diffusion:
% 2D coronal
if ~isempty(answer{2})
    ser = C(str2double(strsplit(answer{2})),1);
    [I,~,fov,p] = readBrukerMRI(fullfile(studydir,ser{1},'pdata','1','2dseq'));
    fov = fov([2,1,3]);
%     I = permute(I,[2,1,3,4]);
    [d(1),d(2),d(3),nvec] = size(I);
    if interpchk
        v = fov./d;
        vd = (d-1)/2 .* v;
        [X,Y,Z] = meshgrid(-vd(2):v(2):vd(2),-vd(1):v(1):vd(1),-vd(3):v(3):vd(3));
        Isave = zeros([d0,nvec]);
        for i = 1:nvec
            Isave(:,:,:,i) = interp3(X,Y,Z,I(:,:,:,i),Xq,Yq,Zq,interpm,0);
        end
        fov = fov0;
    else
        Isave = I;
    end
    adc = zeros(d0);
    for i = 1:3
        adc = adc + log(Isave(:,:,:,1)./Isave(:,:,:,i+1)) / (diff(p.PVM_DwEffBval([1,i+1]))) / 3;
    end
    adc(isnan(adc)) = 0;
    saveMHD(fullfile(svdir,sprintf('%s.mhd',fnout)),...
        flip(cat(4,Isave(:,:,:,1),mean(Isave(:,:,:,2:end),4),adc),1),...
        strcat('DWI_',{'T2w','HighB','ADC'}),fov);
end

%% MT
% assume acquired in 3D coronal
if ~isempty(answer{3})
    ser = C(str2double(strsplit(answer{3})),1);
    MToff = C(str2double(strsplit(answer{3})),10);
    MTi = strcmp(MToff,'off');
    [I,~,fov] = readBrukerMRI(fullfile(studydir,ser{1},'pdata','1','2dseq'));
    [I(:,:,:,2),~] = readBrukerMRI(fullfile(studydir,ser{2},'pdata','1','2dseq'));
    fov = fov([2,1,3]);
    [d(1),d(2),d(3),nvec] = size(I);
    if interpchk
        v = fov./d;
        vd = (d-1)/2 .* v;
        [X,Y,Z] = meshgrid(-vd(2):v(2):vd(2),-vd(1):v(1):vd(1),-vd(3):v(3):vd(3));
        Isave = zeros([d0,nvec]);
        for i = 1:nvec
            Isave(:,:,:,i) = interp3(X,Y,Z,I(:,:,:,i),Xq,Yq,Zq,interpm,0);
        end
        fov = fov0;
    else
        Isave = I;
    end
    MTR = Isave(:,:,:,~MTi)./Isave(:,:,:,MTi);
    MTR(Isave(:,:,:,MTi)<(NoiseLevel(Isave(:,:,round(d(3)/2),MTi))*20)) = 0; % Mask to SNR>20
    saveMHD(fullfile(svdir,sprintf('%s.mhd',fnout)),...
        flip(flip(cat(4,Isave(:,:,:,MTi),Isave(:,:,:,~MTi),MTR),3),1),...
        strcat('MT_',{sprintf('%i',MToff{~MTi}),'off','MTR'}),fov);
end



