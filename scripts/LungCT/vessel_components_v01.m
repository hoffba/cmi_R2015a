function tempResults = vessel_components_v01;

selpath = uigetdir(pwd,'Select Home directory for 5yr Reg Data');
dir_path = dir(selpath);
clear prmFile data_path ID_names

dir_path(contains({dir_path.name},{'output'})) = [];
dir_path(contains({dir_path.name},{'.'})) = [];

NumCases = length(dir_path);

tStart = tic;

for i=1:NumCases
    ID_files = dir(fullfile(dir_path(i).folder,dir_path(i).name)); % list of all imaging data in case folder
    ID_names = dir_path(i).name; % case folder name (e.g. 10000G)
    disp(ID_names)
    
    disp('Load yr0 tPRM and label and yr5 CT scans')
    % Data path
    data_path = fullfile(fullfile(dir_path(i).folder,dir_path(i).name));
    dir_files = {ID_files.name}; % create list of cells for files
    
    %% Load baseline Exp_label
    insSegfile = dir_files(contains(dir_files,'INSP')&contains(dir_files,'COPD.')...
        &~contains(dir_files,'warped')&contains(dir_files,'segmentation')...
        &~contains(dir_files,'.mat'));
    
    insSeg = single(niftiread((fullfile(data_path,insSegfile{:}))));
    
    lobes = unique(insSeg); lobes(lobes==0)=[]; % set ExpL00 as vector "Lobes"
    
    %%
    Vessfile = dir_files(contains(dir_files,'binVessel')&~contains(dir_files,'INS')...
        &~contains(dir_files,'COPD2.')&~contains(dir_files,'warped')&~contains(dir_files,'segmentation')...
        &~contains(dir_files,'.mat'));
    
    vesselSeg = niftiread((fullfile(data_path,Vessfile{:})));
    
    %% Vvessel Search for *.mat file if it does not exist cacluate and save
    
    Vvesselfile = dir_files(contains(dir_files,'.mat'));
    
    if ~isempty(Vvesselfile)
        'Load PreSaved Vvessel'
        Vves = load((fullfile(data_path,Vvesselfile{:})));
        
        Vvessel = Vves.x.img;
        voxsz = Vves.x.voxsz;
    else
        
        disp('Need to load Exp05R and Ins05R to generate PRM then Vvessel')
        % PRM and Vvessel at yr 5
        % 5yr ExpR
        
        % Ins for voxsz
        Insfile = dir_files(contains(dir_files,'INS')&contains(dir_files,'COPD.ct')...
            &~contains(dir_files,'warped')&~contains(dir_files,'segmentation')...
            &~contains(dir_files,'.mat'));
        
        Ins_info = niftiinfo((fullfile(data_path,Insfile{:})));
        voxsz = Ins_info.PixelDimensions;
        
        %%    % Generate PRM for baseline
        
        disp('Generate Vvessel')
        mask = logical(insSeg);
        
        % Topology Analysis (Hoff et al., 2017 Scientific Reports)
        
        str_MF = {'V_e4','S_e4','B_e5','X_e5'};
        
        MFi_final = 1; % this is for volume density
        jj_final = 1; % this is for PRM class (default emph)
        
        
        for jj = 1:jj_final
            MF_mask = ismember(vesselSeg,1);
            
            p = minkowskiFun(logical(MF_mask),'thresh',1,'tmode','==','n',[10,10,10],'gridsp',[5,5,5],'voxsz',voxsz,'mask',logical(mask));
            Vvessel = zeros(size(MF_mask));
            
            mask = []; MF_mask = [];
            
            for MFi = 1:MFi_final
                scale=[1e4 1e4 1e5 1e5];
                
                Vvessel(:,:,:,MFi) = (grid2img(p.MF(1,MFi,:),p.ind,p.mask,3,1).*scale(MFi));
            end
        end
        
        VvesselFinal.img = single(Vvessel);
        VvesselFinal.voxsz = voxsz;
        
%         assignin('base','Vves',VvesselFinal);
        'Save Vvessel'
        parsave(fullfile(data_path, cell2mat([ID_names,'_Vvessel_',str_MF(MFi)])),VvesselFinal)
        %         save(fullfile(data_path, cell2mat([ID_names,'_COPD2_',str_PRMclass,'_',str_MF(MFi)])),'VvesselYr5');
    end
    
    %% Determine components
    'Component Analysis'
    skel=vesselSeg>0;
    CC=bwconncomp(skel);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [L,I]=sort(numPixels,'descend'); % L is the number of elements within a component
    Pixlist=CC.PixelIdxList(I);
    
    B = 100*cumsum(L)/sum(L);
    
    finalPix = Pixlist(B<=95); % find cumsum of components <95% of total vessel numel.
    
    num_ep = [];
    num_vox = [];
    lobeNum = [];
    meanRes = [];
    numComp = [];
    numEp = [];
    VDvess = [];
    
    for i = 1:length(finalPix)
        
        BW_body2 = single(zeros(size(vesselSeg)));
        BW_body2(finalPix{i}) = 1;
        
        A_cl=bwskel(logical(BW_body2));
        uniqLobTemp = unique(BW_body2.*insSeg);
        % A_cl_vox=im2coords3Dbin(A_cl); % Alex code
        
        [~,node,link]=Skel2Graph3D(A_cl,0); % filesharing matlab
        %     figure(1000);plot_graph(node,link); % Alex code
        
        tbl = struct2table(node);
        num_ep = cat(1,num_ep,sum(tbl.ep));
        lobeNum = cat(1,lobeNum, uniqLobTemp(2)); % lobe number
        
    end
    
    numVoxelsPerComp = L(B<=95)';
    tempResults = table(lobeNum, numVoxelsPerComp, num_ep);
    for i = 1:5 % number of lobes
        if ~isempty(numel(tempResults{tempResults{:,1}==i,1}))
            meanRes = cat(1,meanRes, mean(tempResults{tempResults{:,1}==i,1:2},1));
            numEp = cat(1,numEp, sum(tempResults{tempResults{:,1}==i,3},1));
            numComp = cat(1,numComp, numel(tempResults{tempResults{:,1}==i,1}));
            VDvess = cat(1, VDvess, mean(Vvessel(insSeg==i),'all'));
        else
            meanRes = cat(1, meanRes, [0 0]);
            numEp = cat(1, numEp, 0);
            numComp = cat(1, numComp, 0);
            VDves = cat(1, VDves, 0);
        end
    end
    
    
    %% final table config.
    finalResults = splitvars(table(meanRes, numEp, numComp, VDvess));
    
    finalResults.Properties.VariableNames = [{'Lobe'},{'Number Voxels (<95%)'}, ...
        {'Number EndPoints (<95%)'}, {'Number of Components (<95%)'}, {'Vessel Volume Density (Total)'}];
    
    if ~exist(fullfile(selpath,'outputFinal'),'dir')
        mkdir(fullfile(selpath,'outputFinal'))
        mkdir(fullfile(selpath,'outputComponents'))
    end
    parsave(fullfile(selpath,'outputFinal',['Data',ID_names,'.mat']),finalResults);
    parsave(fullfile(selpath,'outputComponents',['Data',ID_names,'components.mat']),tempResults);
end

tEnd = toc(tStart);
formatSpec = 'Total time : %4.1f min\n';
fprintf(formatSpec,tEnd/60)