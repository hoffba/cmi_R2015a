function fwi_Bruker(studydir)
% Assumed acquisition:



interpchk = true;
interpm = 'cubic';

[~,studystr] = fileparts(studydir);
studystr = strsplit(studystr,'_');
fnout = [strjoin(studystr(3:end-2),'-'),'_',studystr{1}];
svdir = fullfile(studydir,studystr{1});
if ~isdir(svdir)
    mkdir(svdir);
end

% Catalog study directory:
[C,hfig] = mriCatalog(studydir,1);

% Find relevant images:
flag = [1,1,1];
ind = {find( cellfun(@(x,y)~isempty(strfind(x,'_Fat'))&&~strcmp(y,'FieldMap'),C(:,2),C(:,3)) )',...
        find(cellfun(@(x)~isempty(strfind(x,'_DtiStandard')),C(:,2)))',...
        find(cellfun(@(x)~isempty(strfind(x,'_MT')),C(:,2)))'};
% answer = inputdlg({'Fat:','Diffusion:','MT'},'MF processing',1,cellfun(@num2str,ind,'UniformOutput',false));
answer = cellfun(@num2str,ind,'UniformOutput',false);
delete(hfig);

%% Fat:
% 3D coronal
if ~isempty(answer{1})
    ser = C(str2double(strsplit(answer{1})),1);
    if ~flag(1)
        % Only load image for geometry:
        [timg,~,fov,p] = readBrukerMRI(fullfile(studydir,ser{1},'fid'));
        [d(1),d(2),d(3),~] = size(timg);
        info.dircos = p.SliceOrient;
        info.slcpos = p.SlicePos;
    else
        for i = 1:length(ser)
            fprintf('Loading series: %s\n',ser{i});
            [timg,~,tfov,p] = readBrukerMRI(fullfile(studydir,ser{i},'fid'));
            [td(1),td(2),td(3),~] = size(timg);
            if i==1
                I = coilCombine(timg);
                fov = tfov;
                te = p.EffectiveTE;
                f0 = p.PVM_FrqRef(1);
                d = td;
                info.dircos = p.SliceOrient;
                info.slcpos = p.SlicePos;
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
                        'SUBSAMPLE',6,...
                        'DO_OT',1,...
                        'LMAP_POWER',2,...
                        'lambda',0.05,...
                        'LMAP_EXTRA',0.05,...
                        'TRY_PERIODIC_RESIDUAL',0);
        % 2D slice by slice:
        di = 1;
        dord = [setdiff(1:3,di),di];
        I = permute(I,[dord,4,5]);
        F = zeros(d(dord)); W = F; R2s = F; fmap = F;
        for i = 1:d(di)
            fprintf('  Slice %u / %u\n',i,d(di));
            tdata.images = I(:,:,i,:,:);
            results = fw_i2cm1i_3pluspoint_hernando_graphcut(tdata,params);
            % Collect results: (F/W are swapped for some reason)
            W(:,:,i) = results.species(1).amps;
            F(:,:,i) = results.species(2).amps;
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
            cat(4,FF,abs(W),abs(F),R2s,fmap,mean(abs(I),5)),...
            strcat('FWI_',{'FatPct','Water','Fat','R2star','FieldMap','MGEmean'}),fov,info);
        clear F W R2s fmap FF I
    end
else
    error('MGE_Fat images not found.');
end

%% Determine interpolation image geometry:
% Assumes all images are not oblique
if interpchk
    fov0 = fov;
    d0 = d;
    [Xq,Yq,Zq] = getImgCoords(d0,fov0./d0,info.dircos,info.slcpos);
end

%% Diffusion:
% 2D coronal
if ~isempty(answer{2}) && flag(2)
    ser = C(str2double(strsplit(answer{2})),1);
    fprintf('Processing Diffusion series: %s\n',ser{1});
    [I,~,fov,p] = readBrukerMRI(fullfile(studydir,ser{1},'pdata','1','2dseq'));
    [d(1),d(2),d(3),nvec] = size(I);
    if interpchk
        [X,Y,Z] = getImgCoords(d,fov./d,p.SliceOrient,p.SlicePos);
        F = scatteredInterpolant(X(:),Y(:),Z(:),reshape(I(:,:,:,1),[],1),'linear','nearest');
        Isave = zeros([d0,nvec]);
        for i = 1:nvec
            F.Values = reshape(I(:,:,:,i),[],1);
            Isave(:,:,:,i) = F(Xq,Yq,Zq);
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
        cat(4,Isave(:,:,:,1),mean(Isave(:,:,:,2:end),4),adc*1000),...
        strcat('DWI_',{'T2w','HighB','ADC'}),fov,info);
end

%% MT
% assume acquired in 3D coronal
if ~isempty(answer{3}) && flag(3)
    ser = C(str2double(strsplit(answer{3})),1);
    fprintf('Processing MT series: %s and %s\n',ser{:});
    MToff = C(str2double(strsplit(answer{3})),10);
    MTi = strcmp(MToff,'off');
    [I,~,fov,p] = readBrukerMRI(fullfile(studydir,ser{1},'pdata','1','2dseq'));
    [I(:,:,:,2),~] = readBrukerMRI(fullfile(studydir,ser{2},'pdata','1','2dseq'));
    [d(1),d(2),d(3),nvec] = size(I);
    if interpchk
        [X,Y,Z] = getImgCoords(d,fov./d,p.SliceOrient,p.SlicePos);
        F = scatteredInterpolant(X(:),Y(:),Z(:),reshape(I(:,:,:,1),[],1),'linear','nearest');
        Isave = zeros([d0,nvec]);
        for i = 1:nvec
            F.Values = reshape(I(:,:,:,i),[],1);
            Isave(:,:,:,i) = F(Xq,Yq,Zq);
        end
        fov = fov0;
    else
        Isave = I;
    end
    MTR = ( Isave(:,:,:,MTi) - Isave(:,:,:,~MTi) ) ./ Isave(:,:,:,MTi);
    MTR(Isave(:,:,:,MTi)<(NoiseLevel(Isave(:,:,round(d(3)/2),MTi))*20)) = 0; % Mask to SNR>20
    saveMHD(fullfile(svdir,sprintf('%s.mhd',fnout)),...
        flip(cat(4,Isave(:,:,:,MTi),Isave(:,:,:,~MTi),MTR),2),...
        strcat('MT_',{sprintf('%i',MToff{~MTi}),'off','MTR'}),fov,info);
end



