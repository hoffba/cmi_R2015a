% CMI script
function runCTlongitudinalPRM(cmiObj)
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



%% Begin find PRM image sets and process
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
[fnames] = parse_study_to_PRM_fnamesguilevels();
assignin('base','fnames',fnames);  % for debugging only


%% Compute prm for filenames loaded

if isempty(fnames)
    fprintf('ERROR: no filenames to process from excel file. Exiting script now.\n');
    return;
end

%% Initialize the analysis parameters - for now, leave at default cmi settings


%% Loop around all selected data sets collecting prm data
for i = 1:length(fnames)
    
    % Feedback for user
    fout = fnames{i}; 
    fout = fout{1};
    [fdir, fout] = fileparts(fout);

    
    % Note: file names are passes as cell array with FULL PATH in each
    fprintf('Begin: %s\n',fout); tic
    
    cmiObj.loadImg(0,fnames{i}(1:2));
    cmiObj.loadMask(fnames{i}{3}); % current cmi program expects character input, not cell,and only one mask   
    cmiObj.img.imgScale([1,2],[1,1],[-1024,-1024]);
    
    fprintf('Done loading %s: ',fout); 
    toc

    
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
    
    % Accumulate stats
    prm_excel(i,:) = outvals;
    
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
fprintf(fid,'Patient,');
for ii=1:length(outlabels)
    fprintf(fid,'%s,',outlabels{ii});
end; fprintf(fid,'\n');

for tt=1:length(fnames)
    sout = fnames{tt}; sout = sout{1};
    [~, sout] = fileparts(sout);
    fprintf(fid,'%s,',sout);
    for ii=1:size(prm_excel,2)
        fprintf(fid,'%6.2f,',prm_excel(tt,ii)); 
    end; fprintf(fid,'\n'); 
end % filenames
fclose(fid);

% Save prm values for later work
assignin('base','out_vals',outvals);

% Feedback for user
fprintf('PRM: summary .csv file %s for all cases saved\n',out_tname);


%% Done

disp('Done with script batch_PRM_excel_cmi!')

