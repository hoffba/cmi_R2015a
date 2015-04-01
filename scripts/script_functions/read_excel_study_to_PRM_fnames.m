function [ fnames, exp_mean_blood, exp_mean_air, ins_mean_blood, ins_mean_air ] = read_excel_study_to_PRM_fnames(  )
%read_excel_study_to_PRM_fnames Make list of filenames from spreadsheet.
% Input: 
%    Interactive: 
%         name of Excel file of study information
% Output:
%    fnames - rows of cells with <Exp> and <InsR> and <ExpVOI>
%    exp_mean_blood, exp_mean_air, ins_mean_blood, ins_mean_air
%
% Description:
% Select excel file for study and read in patient/timepoints/mean blood & air HU
% Each Line:
%     Patient_Id Timept Exp_Mean_blood Exp_Mean_air Ins_Mean_blood Ins_Mean_air
% 
% Revision History:
% 12 July 2013 JLB created
% 22 Oct 2013  JLB revised to work with cmi program script for processing
%                  longitudinal data. Reads dir, fnames, and corr factors
%                   all from excell sheet
% 10 Feb 2014   JLB add in error check for correction factors & default


%% Select excel file for study and read in patient/timepoints/mean blood & air HU
%
% Patient_Id Timept Exp_Mean_blood Exp_Mean_air Ins_Mean_blood Ins_Mean_air
% Note: txt will contain columnn headings

[filename,pathname,filterindex] = uigetfile('*.xlsx','Excel study summary file');
if (filterindex ~= 0)
    cd(pathname);
    [num, txt, raw] = xlsread(filename);
else
    disp('Error: no study file selected. Abort');
    return;
end


% In case I want to expand the script to correct
% Add in default for air/blood HU in case not included
if size(num,2)>=4
    exp_mean_blood = num(:,2);
    exp_mean_air = num(:,1);
    ins_mean_blood = num(:,4);
    ins_mean_air = num(:,3);
else
    exp_mean_blood = ones(size(txt,1)-1).*37;
    exp_mean_air =ones(size(txt,1)-1).*-995;
    ins_mean_blood = ones(size(txt,1)-1).*37;
    ins_mean_air = ones(size(txt,1)-1).*-995;    
    fprintf('WARNING: no air or blood HU values read. Defaulting to 37,-995\n');
end


txt = txt(2:end,:); % first row is label, note these are cells

for tt=1:size(txt,1)
 
    fnames{tt}{1} = fullfile(txt{tt,1}, strcat(txt{tt,2},'.fld'));    
    fnames{tt}{2} = fullfile(txt{tt,1}, strcat(txt{tt,3},'.fld'));
    fnames{tt}{3} = fullfile(txt{tt,1}, strcat(txt{tt,4},'.fld'));
    
end % end of for loop with tt for each timepoint

end % function

