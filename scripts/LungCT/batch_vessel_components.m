function batch_vessel_components(studypath)
% Input: casepath = <cellstr> subject folders
%           folders contain: INSP segmentation,

dir_path = dir(studypath);

dir_path = {dir_path.name};
dir_path(contains(dir_path,{'output'})) = [];
dir_path(contains(dir_path,{'.'})) = [];

NumCases = numel(dir_path);

tStart = tic;

parpool(10);
parfor i = 1:NumCases
    ID = dir_path{i}; % case folder name (e.g. 10000G)
    try
        %% Find files in subject folder:
        data_path = fullfile(fullfile(studypath,ID));
        disp(ID)
        
        %% Find necessary files:
        fn_ins = dir(fullfile(data_path,'re_*_INSP_*.ct.nii*'));
        fn_insSeg = dir(fullfile(data_path,'re_*_INSP_*.lobe_segmentation.nii*'));
        fn_vesselSeg = fullfile(data_path,sprintf('%s_binVessel.nii.gz',ID));
        fn_Vvessel = fullfile(data_path,sprintf('%s_binVessel.tprm.Ve4.nii.gz',ID));
        if ~isempty(fn_insSeg) && ~isempty(fn_ins) && exist(fn_vesselSeg,'file')
            %% Load files:
            fprintf('Load INSP segmentation\n');
            insSeg = single(niftiread(fullfile(data_path,fn_insSeg.name)));
            
            disp('Load info from INSP image');
            info = niftiinfo(fullfile(data_path,fn_ins.name));
            voxsz = info.PixelDimensions;
            
            fprintf('Loading binVessel\n');
            vesselSeg = logical(niftiread(fn_vesselSeg));
            
            %% Load/Generate Vvessel
            if exist(fn_Vvessel,'file')
                disp('Loading Vvessel map from file')
                Vvessel = niftiread(fn_Vvessel);
            else
                disp('Calculating Vvessel map...');
                % Topology Analysis (Hoff et al., 2017 Scientific Reports)
                p = minkowskiFun(vesselSeg,'n',[10,10,10],'gridsp',[5,5,5],'voxsz',voxsz,'mask',logical(insSeg));
                Vvessel = grid2img(p.MF(1,1,:),p.ind,p.mask,3,1).*1e4;
                niftiwrite(int16(Vvessel),fn_Vvessel(1:end-3),info,'Compressed',true);
            end
            
            %% Vessel Component Analysis
            disp('Component Analysis')
            CC = bwconncomp(vesselSeg);
            numPixels = cellfun(@numel,CC.PixelIdxList);
            [L,I] = sort(numPixels,'descend'); % L is the number of elements within a component
            Pixlist = CC.PixelIdxList(I);
            
            B = 100*cumsum(L)/sum(L);
            
            finalPix = Pixlist(B<=95); % find cumsum of components <95% of total vessel numel.
            
            nPix = numel(finalPix);
            lobeNum = zeros(nPix,1);
            num_ep = zeros(nPix,1);
            for ii = 1:nPix
                
                BW_body2 = single(zeros(size(vesselSeg)));
                BW_body2(finalPix{ii}) = 1;
                
                A_cl = bwskel(logical(BW_body2));
                uniqLobTemp = unique(BW_body2.*insSeg);
                % A_cl_vox=im2coords3Dbin(A_cl); % Alex code
                
                [~,node,~] = Skel2Graph3D(A_cl,0); % filesharing matlab
                %     figure(1000);plot_graph(node,link); % Alex code
                
                lobeNum(ii) = uniqLobTemp(2); % lobe number
                if ~isempty(node)
                    num_ep(ii) = sum([node.ep]);
                end
            end
            
            numVoxelsPerComp = L(B<=95)';
            tempResults = table(lobeNum, numVoxelsPerComp, num_ep);
            lobeNum = 1:5;
            meanRes = zeros(5,1); numEp = zeros(5,1); numComp = zeros(5,1); VDvess = zeros(5,1);
            for ii = lobeNum % number of lobes
                ind_lobe = tempResults.lobeNum==ii;
                if ~isempty(numel(tempResults{ind_lobe,1}))
                    meanRes(ii) = mean(tempResults.numVoxelsPerComp(ind_lobe));
                    numEp(ii) = sum(tempResults.num_ep(ind_lobe));
                    numComp(ii) = nnz(ind_lobe);
                    VDvess(ii) = mean(Vvessel(insSeg==ii),'all');
                end
            end
            
            %% final table config.
            finalResults = table(lobeNum', meanRes, numEp, numComp, VDvess);
            
            finalResults.Properties.VariableNames = [{'Lobe'},{'Number Voxels (<95%)'}, ...
                {'Number EndPoints (<95%)'}, {'Number of Components (<95%)'}, {'Vessel Volume Density (Total)'}];
            
            if ~exist(fullfile(studypath,'outputFinal'),'dir')
                mkdir(fullfile(studypath,'outputFinal'))
                mkdir(fullfile(studypath,'outputComponents'))
            end
            writetable(finalResults,fullfile(studypath,'outputFinal',['Data_',ID,'.csv']));
            writetable(tempResults,fullfile(studypath,'outputComponents',['Data_',ID,'_components.csv']));
            %         parsave(fullfile(selpath,'outputFinal',['Data',ID,'.mat']),finalResults);
            %         parsave(fullfile(selpath,'outputComponents',['Data',ID,'components.mat']),tempResults);
        end
    catch err
        fprintf(1,'ID: %s ; %s\n',ID,err.identifier);
        fprintf(1,'There was an error! The message was:\n%s',e.message);
    end
end

tEnd = toc(tStart);
formatSpec = 'Total time : %4.1f min\n';
fprintf(formatSpec,tEnd/60)