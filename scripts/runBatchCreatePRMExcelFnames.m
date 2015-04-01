% CMI script
function runbatchCreatePRMExcelFnames(cmiObj)
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
% 9 Jan 2014   jlb fixed error check for .csv filename
% 10 Jan 2014  jlb begin to add in check option for full filenames & up to 5 files
%                   (fix in subroutine that parses filenames first - done)
%      
% 3 Feb 2014  mb make a check for files past 3 


%% Output file for .csv filename info
[out_tname, out_tpath] = uiputfile('*.csv','Select output filename .csv file');
if isequal(out_tpath,0) || isequal(out_tname,0)
    fprintf('Error: no .csv file name selected. Exiting script.\n');
    return;
end
%% Load  names plus additional data from an Excel file
[fnames pathFlag] = parse_study_to_PRM_fnamesguilevels5();
if isempty(fnames)
    fprintf('ERROR: found no filenames to process. Exiting script now.\n');
    return;
end
assignin('base','fnames',fnames);  % for debugging only


%% Output all PRM fnames manually to .csv file


% Open file
csvfname= fullfile(out_tpath,out_tname);
fid = fopen(csvfname,'w');
if fid==-1
    fprintf('****ERROR opening %s\n',out_tname);
    return;
end

% output ['Patient' outlabels] and [patient_id outvals]
fprintf(fid,'PatientBaseDir,');
for m=1:5
    fprintf(fid,'File%d,',m);
end
fprintf(fid,'\n');

for tt=1:length(fnames)
    sout = fnames{tt}; 
    for m=1:length(sout)
        if isempty(sout{m})
            fprintf(fid,'.,')
        else 
            [basedir, sout1] = fileparts(sout{m});
            if m==1
                fprintf(fid,'%s,',basedir);
            end
            if pathFlag
                fprintf(fid,'%s,',sout{m});
            else
                fprintf(fid,'%s,',sout1);
            end
        end
    end
    fprintf(fid,'\n');
end % filenames
fclose(fid);

% Save prm values for later work
%assignin('base','out_vals',outvals);

% Feedback for user
fprintf('PRM: summary .csv file %s for all cases saved\n',out_tname);


%% Done

disp('Done with script runBatchCreatePRMExcelFnames!')


