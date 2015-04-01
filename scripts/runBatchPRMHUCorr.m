% CMI script
function runBatchPRMHUCorr(cmiObj)
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
% 7 Feb 2014   jlb add option to process dicom files and not subtract 1024
% 10 Feb 2014  jlb add error check for correction factors and default
%                  (note - mostly in function read...)
%                  add write of full 3 filenames prm generated from
% 3 Mar 2014   jlb add in option to erode VOI before processing
            
%% Set PRM options & ERODE options
disp('Please set the PRM options before you start!');
cmiObj.setPRMopts();

disp('Please set the ERODE flag for VOI exploring');
button = questdlg('ERODE VOI for experiment?');
switch button
    case 'Yes'
        erode_flag = 1;
        
    case 'No'
        erode_flag = 0;
        
    case 'Cancel'
        erode_flag = 0;  
end


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
% Note: Add in check for reading full filenames or 
% patient directories/single filenames
button = questdlg('Do filenames contain full directory paths?');
switch button
    case 'Yes'
        % Full filename order: <ignored> <VOI> <Exp> <Ins> (<E> <I>)
        % 
        [fnames, exp_mean_blood, exp_mean_air, ins_mean_blood, ins_mean_air] = read_excel_study_to_PRM_fullfnames();
        % files in a different order, extract the ones we need
        for i=1:length(fnames)
            fnamesnew{i}{1} = fnames{i}{3};
            fnamesnew{i}{2} = fnames{i}{4};
            fnamesnew{i}{3} = fnames{i}{2};
        end
        clear fnames;
        fnames = fnamesnew;     
        
    case 'No'
        [fnames, exp_mean_blood, exp_mean_air, ins_mean_blood, ins_mean_air] = read_excel_study_to_PRM_fnames();
        
    case 'Cancel'
        [fnames, exp_mean_blood, exp_mean_air, ins_mean_blood, ins_mean_air] = read_excel_study_to_PRM_fnames();      
end
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
    offset = [-1024 -1024];
    try 
        dicominfo(fnames{i}{1});
        dflag1 = 1;
        offset(1) = 0;
    catch
        dflag1 = 0;
    end
    try 
        dicominfo(fnames{i}{2});
        dflag2 = 1;
        offset(2) = 0;
    catch
        dflag2 = 0;
    end    
    stat = cmiObj.loadImg(false,fnames{i}(1:2));
    if ~stat || cmiObj.img.dims(4)<2
         % Error - didn't load 2 same size images
        % set prm to zeros and break to next set
        fprintf('Error loading at least one of file %s: moving to next set',fout);
        if i~=1
            outvals = outvals .* 0;
            prm_excel(i,:) =  outvals;
        end
        
    else
        
        cmiObj.loadMask(fnames{i}{3}); % current cmi program expects character input, not cell,and only one mask
        cmiObj.img.imgScale([1,2],[1,1],offset);
        
        % Erode VOIs if you are doing this experiment
        if erode_flag
            cmiObj.img.mask.morph('erode',[1 1]); % erode by one in xy and one in z
        end
        
        % Files loaded. Compute PRM
        
        fprintf('Done loading %s: ',fout);
        toc
        
        % HU correction
        fprintf('Start HU correcting');
        aa = [exp_mean_air(i) ins_mean_air(i)];
        bb = [exp_mean_blood(i) ins_mean_blood(i)];
        vec = [1 2];
        cmiObj.img.imgHUCorrect(vec,aa,bb);
        %refAirBlood = [-1000 50];
        %cmiObj.img.imgHUCorrect(vec,aa,bb,refAirBlood);        
        
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
        
        % Combine prmvals with stats
        outlabels = [prmlabels statlabels];
        outvals = [prmvals statvals];
        
        % Feedback for user
        for ii=1:length(outlabels)
            fprintf('%s ',outlabels{ii});
        end; fprintf('\n');
        for ii=1:length(outvals)
            fprintf('%6.2f ',outvals(ii));
        end; fprintf('\n');
        
    end % if 2 files to process
    
    % Accumulate stats
    prm_excel(i,:) = outvals;
    
    %Print off .mat file for each patient - don't save .mat files for now
    % JB 10 feb 2014
    %[pdir, pname, pext] = fileparts(fnames{i}{2});
    %prm_mat= cmiObj.img.prm.mat;
    %save(fullfile(pdir, [pname '_prm.mat']),'prm_mat');
    
    
    % create PRM subdirctory full path and name
    %pfdir = fullfile(fdir,prm_dir);
    
end % for each triple of files



%% Output all PRM values manually to .csv file


% Open file
csvfname= fullfile(out_tpath,out_tname);
fid = fopen(csvfname,'w');
if fid==-1
    fprintf('****ERROR opening %s\n',out_tname);
    return;
end

% output ['Patient' outlabels] and [patient_id outvals]
fprintf(fid,'File_Exp,File_Ins,File_VOI,');
for ii=1:length(outlabels)
    fprintf(fid,'%s,',outlabels{ii});
end; 
fprintf(fid,'\n');

for tt=1:length(fnames)
    sout = fnames{tt}; sout = sout{1};
    [pdir, sout] = fileparts(sout);
    fprintf(fid,'%s,%s,%s,',fnames{tt}{1},fnames{tt}{2},fnames{tt}{3});
    for ii=1:size(prm_excel,2)
        fprintf(fid,'%6.2f,',prm_excel(tt,ii));
    end; 
    fprintf(fid,'\n');
end % filenames
fclose(fid);

% Save prm values for later work
assignin('base','out_vals',outvals);

% Feedback for user
fprintf('PRM: summary .csv file %s for all cases saved\n',out_tname);


%% Done

disp('Done with script batch_PRM_excel_cmi!')


