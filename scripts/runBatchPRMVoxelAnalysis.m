% CMI script
function runBatchPRMVoxelAnalysis(cmiObj)
%   Input: cmiObj = CMIclass object containing current settings
%   Input: cmiObj = CMIclass object containing current settings
%   Input: cmiObj = CMIclass object containing current settings
%
% Performs a batch analysis of Lung-CT PRM data
% 1. Read in excel file and generate name triplets
%     - all must be in the same folder
%     - tested on FLDs
% 2. Cycles through each set and calculates PRM statistics
%       a) calculate 6 color PRM to have all the values accessible
%           and save to array
%       b) calculate emphy, gas trapping volumes, lung volumes?
%       c) save .csv file of all PRMs calculated
%
% Revision history:
% 12 July 2013 jlb converted from .mat/user cycle to reading Excel file
% 22 Sep  2013 jlb converted to cmi program launched script that uses all 
%                  cmi functionality for doing computations
% 8 Oct 2013   jlb added in option for specifying file formatting
% 13 Jan 2014  jlb transformed from batch PRM to batch voxel-wise analysis
            
%% Set PRM options
disp('Please set the PRM options before you start!');
cmiObj.setPRMopts();


%% Begin find PRM image sets and process
HU_CORRECT = 1;
MAX_PRM_STATS = 15;
EMPHYSEMA_THRESHOLD = -950; % inspiration or image 2
GAS_TRAPPING_THRESHOLD = -856; % expiration or image 1
prm_dir = 'PRM'; % Set this to the subdirectory for saving PRM results for each case

%% Output file for .csv prm info
[out_tname, out_tpath] = uiputfile('*.csv','Select output .csv file');
if isequal(out_tpath,0) || isequal(out_tname,0)
    fprintf('Error: no .csv file name selected. Exiting script.\n');
    return;
end
%% Load  names plus additional data from an Excel file
[fnames, exp_mean_blood, exp_mean_air, ins_mean_blood, ins_mean_air] = read_excel_study_to_PRM_fullfnames();
assignin('base','fnames',fnames);  % for debugging only


%% Compute prm for filenames loaded

if isempty(fnames)
    fprintf('ERROR: no filenames to process from excel file. Exiting script now.\n');
    return;
end

%% Loop around all selected data sets collecting prm data
for i = 1:length(fnames)
    
    % Feedback for user
    [fdir, fout] = fileparts(fnames{i}{1});
    fprintf('Begin: %s\n',fout); tic
    
    for m=1:2
        % Load images
        % 1:Patient Dir 2:VOI 3-6:Exp InsR ExpT2RExpT1 InsRRExpT1
        if m==1 
            cmiObj.loadImg(0,fnames{i}(3:4)); 
            alloutvals = [];
        else
            cmiObj.loadImg(0,fnames{i}(5:6)); 
        end;
        cmiObj.loadMask(fnames{i}{2}); % current cmi program expects character input, not cell,and only one mask
        cmiObj.img.imgScale([1,2],[1,1],[-1024,-1024]);
        
        fprintf('Done loading %s: ',fout);
        toc
        
        % HU correct images
        fprintf('Start HU correcting...\n');
        aa = [exp_mean_air(i) ins_mean_air(i)];
        bb = [exp_mean_blood(i) ins_mean_blood(i)];
        vec = [1 2];
        cmiObj.img.imgHUCorrect(vec,aa,bb);
        
        % Calculate PRM with:
        %       Expiration as x-axis (vec 1 is default)
        %       Inspiration as y-axis (vec 2 passed in)
        fprintf('Start computing PRM...\n'); tic
        
        cmiObj.img.calcPRM(2);
        [prmlabels prmvals] = cmiObj.img.prm.getStats();
        fprintf('Done with PRM: '); toc
        
        % Find additional stats
        voxvol = prod(cmiObj.img.voxsz);
        npm = numel(find(cmiObj.img.mask.mat));
        np_prm = numel(find(cmiObj.img.prm.mat));
        
        [vMean,vStD,vMed,vVol,nVox] = cmiObj.img.getStats; % This is working now
        
        statlabels(1) = {'VOI_Volume'};
        statvals(1) = vVol(1);
        
        statlabels(2) = {'Emphy_Index_pseudo'};
        statvals(2) = 100*numel(find(cmiObj.img.mask.mat & cmiObj.img.mat(:,:,:,2)<=-950))/npm;
        
        statlabels(3) = {'Gas_Trapping_Index_pseudo'};
        statvals(3) = 100*numel(find(cmiObj.img.mask.mat & cmiObj.img.mat(:,:,:,1)<=-856))/npm;
        
        statlabels(4) = {'Mean_1'};
        statvals(4) = vMean(1);
        
        statlabels(5) = {'Median_1'};
        statvals(5) = vMed(1);
        
        statlabels(6) = {'Mean_2'};
        statvals(6) = vMean(2);
        
        statlabels(7) = {'Median_2'};
        statvals(7) = vMed(2);
        
        % Combine prmvals with stats - not for now
%         outlabels = [prmlabels statlabels];
%         outvals = [prmvals statvals];
         outlabels = [prmlabels];
         outvals = [prmvals];
        
        % Feedback for user
        for ii=1:length(outlabels)
            fprintf('%s ',outlabels{ii});
        end; fprintf('\n');
        for ii=1:length(outvals)
            fprintf('%6.2f ',outvals(ii));
        end; fprintf('\n');
        
        % Accumulate stats
        alloutvals = [alloutvals outvals];
        if m==1
            prm1 = cmiObj.img.prm.mat;
        else
            prm2 = cmiObj.img.prm.mat;
        end
    
    end
    prm_excel(i,:) = alloutvals;
    
    % Compute voxel transitions from prm1 and prm2
    prmTransTo = zeros(cmiObj.img.prm.nprm,cmiObj.img.prm.nprm);
    prmTransFrom = zeros(cmiObj.img.prm.nprm,cmiObj.img.prm.nprm); 

    for pp=1:cmiObj.img.prm.nprm
        for qq=1:cmiObj.img.prm.nprm % cycle thru columns of prm
            prmTransTo(pp,qq) = nnz(prm1==qq & prm2==pp);
            prmTransFrom(pp,qq) = nnz(prm1==pp & prm2==qq);
        end
    end
    
    kindex = 1;
    prm4Kappa = zeros(cmiObj.img.prm.nprm*cmiObj.img.prm.nprm,3);
    for pp=1:cmiObj.img.prm.nprm
        for qq=1:cmiObj.img.prm.nprm % cycle thru columns of prm
            prm4Kappa(kindex,1) = pp;
            prm4Kappa(kindex,2) = qq;
            prm4Kappa(kindex,3) = prmTransTo(qq,pp); % cols are from, rows are to
            kindex = kindex+1;
        end
    end   
    prm_excel2(i,:,:) = prmTransTo;
    prm_excel3(i,:,:) = prmTransFrom;
    
    % write out array of label counts for t1 and t2 to .csv file
    writeCSV(fullfile(fdir,sprintf('prm%dKappa.csv',cmiObj.img.prm.nprm)),...
              prm4Kappa,{'PRMLabelT1','PRMLabelT2','VoxelCount'});
    fprintf('Done writing out prmAllLung4Kappa.csv.\n');
    
    
%Print off .mat file for each patient - leave out for now
%     [pdir, pname, pext] = fileparts(fnames{i}{2});
%     prm_mat= cmiObj.img.prm.mat;
%     save(fullfile(pdir, [pname '_prm.mat']),'prm_mat');
    
    

   
end % for each triple of files



%% Output all PRM values manually to .csv file


% Open file
csvfname= fullfile(out_tpath,out_tname);
fid = fopen(csvfname,'w');
if fid==-1
    fprintf('****ERROR opening %s\n',out_tname);
    return;
end

% output ['Patient' outlabels outlables] and [patient_id outvals*2]
fprintf(fid,'Patient,');
for ii=1:length(outlabels)
    fprintf(fid,'PRM_t1_%s,',outlabels{ii});
end; 
for ii=1:length(outlabels)
    fprintf(fid,'PRM_t2_%s,',outlabels{ii});
end; 
fprintf(fid,'\n');

for tt=1:length(fnames)
    fprintf(fid,'%s,',fnames{tt}{1}); % patient base directory here
    for ii=1:size(prm_excel,2)
        fprintf(fid,'%6.2f,',prm_excel(tt,ii)); 
    end; fprintf(fid,'\n'); 
    for ii=1:size(prm_excel2,2)
       fprintf(fid,'%s_prmt1to_prmt2from_%d,',fnames{tt}{1},str2num(cmiObj.img.prm.prmmap{ii,2})); % patient base directory here
       for jj=1:size(prm_excel2,3)
          fprintf(fid,'%d,',prm_excel2(tt,ii,jj)); 
       end
       for jj=1:size(prm_excel3,3)
          fprintf(fid,'%d,',prm_excel3(tt,ii,jj)); 
       end       
       fprintf(fid,'\n');
    end
    
end % filenames


fclose(fid);

% Save prm values for later work
assignin('base','out_vals',outvals);

% Feedback for user
fprintf('PRM: summary .csv file %s for all cases saved\n',out_tname);


%% Done

disp('Done with script batch_PRM_excel_cmi!')



