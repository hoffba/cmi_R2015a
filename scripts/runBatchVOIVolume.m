% CMI script
function runBatchVOIVolume(cmiObj)
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
% 22 Oct 2013  jlb convert to volume only computations (still in progress)


%% Output file for .csv prm info
[out_tname, out_tpath] = uiputfile('*.csv','Select output .csv file');
if isequal(out_tpath,0) || isequal(out_tname,0)
    fprintf('Error: no .csv file name selected. Exiting script.\n');
    return;
end
%% Load  names plus additional data from an Excel file
[fnames] = parse_study_to_VOI_fnamesguilevels();
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
    
    cmiObj.loadImg(0,fnames{i}(1));
    cmiObj.loadMask(fnames{i}{2}); % current cmi program expects character input, not cell,and only one mask
    cmiObj.img.imgScale([1,2],[1,1],[-1024,-1024]);
    
    fprintf('Done loading %s: ',fout);
    toc
    
    
    % Calculate PRM with:
    %       Expiration as x-axis (vec 1 is default)
    %       Inspiration as y-axis (vec 2 passed in)
    fprintf('Start getting stats...\n'); tic
    %[prmlabels prmvals] = cmiObj.img.prm.getStats(); % JB 11-4-2013
    
    % Compute  additional stats
    npm = numel(find(cmiObj.img.mask.mat));
    
    [vMean,vStD,vMed,vVol,nVox] = cmiObj.img.getStats; % This is working now
    
    count = 1;
    statlabels(count) = {'VOI_Volume'};
    if numel(vVol)==0 statvals(count) = -1;
    else statvals(count) = vVol(1); end
    
    count = count+1;
    statlabels(count) = {'Mean_1'};
    if numel(vMean)==0 statvals(count) = -1;
    else statvals(count) = vMean(1); end
    
    count = count + 1;
    statlabels(count) = {'Median_1'};
    if numel(vMed)==0 statvals(count) = -1;
    else statvals(count) = vMed(1); end
    
    
    
    % Combine prmvals with stats
    outlabels = statlabels;
    outvals = statvals;
    
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



%% Output all VOI stat values manually to .csv file


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

