function T = tPRM_aim2_sub(procdir)

tStart = tic;

Vcontour = [0.25 0.5 0.75]*10000; % >25%, 50% and 75%
nc = numel(Vcontour);
edges = reshape(repmat(10*(1:5),nc+2,1) + repmat((0:(nc+1))'-.5,1,5),1,[]);
edges((nc+3):(nc+2):end) = [];

%% Initialize output table
str_MF = ["V_e4","S_e4","B_e5","X_e5"];
str_PRM = ["Norm","fSAD","emph","PD"];
strPRMMF_label = reshape(str_PRM + "_" + str_MF',1,[]);
vars = [{'cellstr','int8','single','int32','single','int32','single','single'},repmat({'single'},1,16) ; ...
    {'ID','Vemph_Contour','Lobe','Lobe_Nvox','Lobe_Vol','Contour_Nvox','Contour_Vol','V_emph_e4_yr5'},strPRMMF_label];
T = table('Size',[nc*5,size(vars,2)],...
          'VariableTypes',vars(1,:),...
          'VariableNames',vars(2,:));
T.Vemph_Contour = repmat(Vcontour',[5,1]);
T.Lobe = reshape(repmat([11 12 13 21 22],[nc,1]),[],1);

fname = dir(procdir); % list of all imaging data in case folder
fname = {fname(3:end).name};
[~,ID] = fileparts(procdir); % case folder name (e.g. 10000G)
T.ID = repmat({ID},15,1);

%% Load baseline Exp_label
fprintf('%s: Load yr0 tPRM and label and yr5 CT scans\n',ID)
fprintf('%s: Load baseline Exp_label\n',ID)
ind = find(~cellfun(@isempty,regexp(fname,'.+_EXP_.+segmentation.nii.gz')),1);
if isempty(ind)
    warning(['Could not find segmentation file to load ... skipping ',ID]);
else
    fn_ExpL00 = fname{ind};
    info = niftiinfo(fullfile(procdir,fn_ExpL00));
    label = single(niftiread(info));
    
    lobes = unique(label); lobes(lobes==0)=[]; % set ExpL00 as vector "Lobes"
    % Relabel lobes so that it is easier to label Vemph threshold in each
    % lobe. e.g. 11 = RUL (10) with Vemph>0.25 (1)
    label(label==11) = 10; % RUL
    label(label==12) = 20; % RLL
    label(label==13) = 30; % RML This is the odd one
    label(label==21) = 40; % LUL
    label(label==22) = 50; % LLL
    lungmask = ismember(label,1:60); % convert label to mask
    
    %% 5yr Vemph Search for *.mat file if it does not exist load Exp05R and Ins05R to generate PRM then Vemph05
    ind = find(~cellfun(@isempty,regexp(fname,'\.mat')),1);
    if ~isempty(ind)
        fprintf('%s: Loading PreSaved Vemph yr 05\n',ID);
        Vemph05 = load(fullfile(procdir,fname{ind}));
        voxsz = Vemph05.voxsz;
        Vemph05 = Vemph05.Vemph05;
    else
        indE = find(~cellfun(@isempty,regexp(fname,'\..+_EXP_.+_COPD2\.warped.nii.gz')),1);
        indI = find(~cellfun(@isempty,regexp(fname,'\..+_INSP_.+_COPD2\.warped.nii.gz')),1);
        if isempty(indE)
            warning('%s: Could not locate EXP image for loading.',ID); return;
        elseif isempty(indI)
            warning('%s: Could not locate INSP image for loading.',ID); return;
        else
            fprintf('%s: Load EXP\n',ID);
            info = niftiinfo(fullfile(procdir,fname{indE}));
            exp = niftiread(info);
            fprintf('%s: Load INSP\n',ID);
            ins = niftiread(fullfile(procdir,fname{indI}));
            voxsz = info.PixelDimensions;
            
            fprintf('%s: Generate yr5 PRM\n',ID)
            exp = medfilt3(double(exp),[3,3,3]);
            ins = medfilt3(double(ins),[3,3,3]);
            prm = PRMclass;
            prm.setModel('AllLung_5Color')
            prm.calcPRM([exp(lungmask),ins(lungmask)],1:2,2,{'Expiration','Inspiration'},lungmask);
            clear exp ins
            prm = prm.mat;
            
            fprintf('%s: Generate yr5 Vemph\n',ID)
            p = minkowskiFun(ismember(prm,[4,5]),'n',[10,10,10],'gridsp',[5,5,5],'voxsz',voxsz,'mask',lungmask,'ord',0);
            clear prm
            Vemph05 = grid2img(p.MF(1,1,:),p.ind,lungmask,3,1) * 1e4;
            save(fullfile(procdir,sprintf('%s_COPD2_emph_V_e4',ID)),'Vemph05','voxsz');
            
        end
    end
    
    %% Generate Vemph contour and calculate lobe and Vemph vol and numel
    fVemph05 = zeros(size(Vemph05));

    % Creat Vemph labels, must be in this order
    for icont = 1:nc
        fVemph05(Vemph05>=Vcontour(icont)) = icont; % V>=0.25
    end
    Vemph05 = Vemph05(lungmask);

    % Unique labels for lobe and contour
    lobeVemph05 = label + fVemph05;

    [voxCont, ~, bin] = histcounts(lobeVemph05(lungmask),edges);
    idxBin = 1:(5*(nc+1));
    % Lobe stats
    T.Lobe_Nvox = reshape(repmat(sum(reshape(voxCont,4,[])),[3,1]),[],1);
    T.Lobe_Vol = T.Lobe_Nvox*prod(voxsz)/1e6;
    % Remove counts/bins outside contours
    ind = 1:(nc+1):numel(voxCont);
    voxCont(ind) = [];
    idxBin(ind) = [];
    % Contour stats
    T.Contour_Nvox = voxCont';
    T.Contour_Vol = voxCont'*prod(voxsz)/1e6;

    % Reshape for easy indexing
    idxBinRe = reshape(idxBin,3,[])';
    
    % Find unique bin numbers
    binList = unique(bin);
    binList(ismember(binList,ind))=[];
    nbin = numel(binList);

    %% load baseline tPRM
    if any(voxCont)
        fprintf('%s: Start working through tPRM baseline\n',ID)
        for iprm = 1:4 % PRM category
            for imf = 1:4 % MF metric
                j = imf + 4*(iprm-1); % index into PRM/MF
                fprintf('%s: Analyzing: %s\n',ID,strPRMMF_label{j});

                % Find tPRM file
                ind = find(~cellfun(@isempty,regexp(fname,['.+_EXP_.+_COPD_',strPRMMF_label{j},'\.nii\.gz'])),1);
                if isempty(ind)
                    warning('%s: tPRM file not found: %s',ID,strPRMMF_label{j});
                else
                    V = zeros(15,1);

                    tPRMData = niftiread(fullfile(procdir,fname{ind}));
                    tPRMData = tPRMData(lungmask);

                    for binI = 1:nbin
                        % Find tPRM yr0 means within region
                        [r,c,~] = find(idxBinRe==binList(binI));
                        [~,loc] = ismember(bin,idxBinRe(r,c:end));
                        ind = idxBin==binList(binI);
                        loc = logical(loc);
                        V(ind) = mean(tPRMData(loc));
                        T.V_emph_e4_yr5(ind) = mean(Vemph05(loc));
                    end

                    % Add results to table:
                    T.(strPRMMF_label{j}) = V;
                end
            end
        end
    else
        fprintf('%s: Contour counts are zero - nothing to analyze\n',ID);
    end
    T = sortrows(T,2);
    
    if ~exist(fullfile(procdir,'output_v04'),'dir')
        mkdir(fullfile(procdir,'output_v04'))
    end
    save(fullfile(procdir,'output_v04',['Data',ID,'.mat']),'T');
        
    fprintf('%s: TOTAL processing time: %s\n',ID,duration(0,0,toc(tStart)));
end