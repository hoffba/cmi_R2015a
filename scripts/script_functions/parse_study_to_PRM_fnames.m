function [ fnames ] = parse_study_to_PRM_fnames(  )
%read_excel_study_to_PRM_fnames Make list of filenames from spreadsheet.
% Input: 
%    Interactive: 
%         name of directory where all patient subdirectories are found
%         name of Excel file of study information
% Output:
%    fnames - rows of cells with <Exp> and <InsR> and <ExpVOI>
%
% Description:
% Select excel file for study and read in patient/timepoints/mean blood & air HU
% Each Line:
%     Patient_Id Timept Exp_Mean_blood Exp_Mean_air Ins_Mean_blood Ins_Mean_air
% 
% Revision History:
% 12 July 2013 JLB created
% 30 Sep  2013 JLB modify to find specific files in directory structure
%                  with error checking


%% Select study directory to process

dataset_folder = uigetdir('/mnt/cmi/projects/CT Lung/clinical/clinical_longitudinal',...
    'Select directory containing all patient subdirectories:');
if (dataset_folder == 0) 
    disp('Error: no study main fld directory specified. Abort')
    return;
end

%% Create subdirectory/file lists from input study summary data

% Filename format
% <>_T#_Exp_bone
% <>_T#_Exp_bone_VOI
% <>_T#_Ins_bone
% <>_T#_Ins_R_bone

savedir = pwd;
ids = dir(dataset_folder);
count = 0;

for i=1:numel(ids)
    % Cycle through IDs in experiment
    if (~strcmp(ids(i).name(1),'.') && ids(i).isdir)
        % process if is directory and not a hidden file
        fprintf('Process ID=%s\n',ids(i).name);
        tpts = dir(fullfile(dataset_folder,ids(i).name));
        
        % cycle through data points for an ID
        for j=1:numel(tpts)
            if (~strcmp(tpts(j).name(1),'.') && tpts(j).isdir)
                % process if is directory and not a hidden file
                fprintf('  data point = %s\n',tpts(j).name);
                cd(fullfile(fullfile(dataset_folder,ids(i).name,tpts(j).name)));
                files = dir('*.fld');
                
                f1 = 0; f2 = 0; f3 = 0;
                for k=1:numel(files)
                    if numel(strfind(files(k).name,'STD_Exp')) &&...
                            ~numel(strfind(files(k).name,'_VOI'))
                        fprintf('       AVS _Ins field file = %s\n',files(k).name)
                        f1 = k;
                    end
                    if numel(strfind(files(k).name,'STD_Ins_R')) &&...
                            ~numel(strfind(files(k).name,'VOI'))
                        fprintf('       AVS _Exp field file = %s\n',files(k).name)
                        f2 = k;
                    end
                    if numel(strfind(files(k).name,'STD_Exp_VOI'))
                        fprintf('       AVS _VOI field file = %s\n',files(k).name)
                        f3 = k;
                    end
                end
                if f1~=0 && f2~=0 && f3~=0
                    count = count + 1;
                    fnames{count}{1} = fullfile(dataset_folder,ids(i).name,tpts(j).name,files(f1).name);
                    fnames{count}{2} = fullfile(dataset_folder,ids(i).name,tpts(j).name,files(f2).name);
                    fnames{count}{3} = fullfile(dataset_folder,ids(i).name,tpts(j).name,files(f3).name);
                end
            end
        end
    end
end

cd(savedir);

% fnames{tt}{1} = fullfile(current_dir, filename_exp);
% fnames{tt}{2} = fullfile(current_dir, filename_insR);
% fnames{tt}{3} = fullfile(current_dir, filename_voi);


end % function

