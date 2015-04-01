function [ fnames pathFlag] = parse_study_to_PRM_fnamesguilevels5(  )
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
% 7  Oct  2013 JLB add in the capability to to change parsinging strings
% 10 Jan  2014 JLB add option for up to 5 filenames & check for full files
% 7  Feb  2014 JLB fix to search all files not just .fld
% 7  Feb  2014 JLB add option for DICOM directory and files
% 12 Feb  2014 JLB debug DICOM error - adding non-DICOM files
% 22 july 2014 JLB add in option for 3 level search

%% Set constant parsing strings

f1format = '_Ins_VOI.fld';
f2format = '_Exp_R.fld';
f3format = '_Ins.fld';
f4format = '~';
f5format = '~';
f6format = '~';
levels = 1;
pathFlag = 0;
dicomFlag = 0;


prompt = {'VOI filename format:' 'First filename format:' 'Second filename format:' 'Third filename format:'...
    'Fourth filename format:' 'Fifth filename format:' 'Directory levels' 'Full filename?' 'DICOM files?'};
title = 'PRM Filename specs';
defAns = {f1format f2format f3format f4format f5format f6format num2str(levels) 'no' 'no'};
options.Resize = 'on';
answer = inputdlg(prompt,title,1,defAns,options);
if isempty(answer)
    disp('Error: cancelling filename specification. Abort')
    fnames = {};
    return;
end

f_format = answer(1:6);
findf = ~strcmp(f_format,'~');
levels = str2double(answer{7});
if isequal(answer{8},'yes')
    pathFlag = 1;
else
    pathFlag = 0;
end
if isequal(answer{9},'yes')
    dicomFlag = 1;
else
    dicomFlag = 0;
end


fprintf('Parse filenames using voi=%s f1=%s f2=%s f3=%s f4=%s f5=%s down %d levels',...
    f_format{1},f_format{2},f_format{3},f_format{4},f_format{5},f_format{6},levels);

%% Select study directory to process

dataset_folder = uigetdir('/mnt/cmi/projects/CT Lung/clinical/clinical_longitudinal',...
    'Select directory containing all patient subdirectories:');
if (dataset_folder == 0)
    disp('Error: no study main fld directory specified. Abort')
    fnames = {};
    return;
end

%% Create subdirectory/file lists from input study summary data


savedir = pwd;
ids = dir(dataset_folder);
count = 0;
fnames = {};

for i=1:numel(ids)
    % Cycle through IDs in experiment
    if (~strcmp(ids(i).name(1),'.') && ids(i).isdir)
        % process if is directory and not a hidden file
        fprintf('Process ID=%s\n',ids(i).name);
        cd(fullfile(dataset_folder,ids(i).name));
        switch levels
            case 1 % data in this directory (for dicoms this will be a dir)
                %files = dir('*.fld');
                files = dir();
                findex = zeros(1,6);
                for k=1:numel(files)
                    for m=1:6
                        if findf(m) &&...
                                numel(strfind(files(k).name,f_format{m}))
                            findex(m) = k;
                        end
                    end
                end
                % print out results of search
                if sum(findex)==0
                    fprintf(' ERROR - 0 files found\n');
                else
                    % Found some or all
                    for m=1:6
                        if findf(m)
                            if findex(m)
                                fprintf('        %s = %s\n',f_format{m},files(findex(m)).name);
                            else
                                fprintf('        %s = ERROR no file found.\n',f_format{m});
                            end
                        end
                    end
                end
                
                if findex(1)>0 && findex(2)>0 && findex(3)>0 % found more than 1 file
                    count = count + 1;
                    if ~dicomFlag
                        fnames{count}{1} = fullfile(dataset_folder,ids(i).name,files(findex(1)).name);
                        fnames{count}{2} = fullfile(dataset_folder,ids(i).name,files(findex(2)).name);
                        fnames{count}{3} = fullfile(dataset_folder,ids(i).name,files(findex(3)).name);
                        if findex(4) fnames{count}{4} = fullfile(dataset_folder,ids(i).name,files(findex(4)).name); end;
                        if findex(5) fnames{count}{5} = fullfile(dataset_folder,ids(i).name,files(findex(5)).name); end;
                        if findex(6) fnames{count}{6} = fullfile(dataset_folder,ids(i).name,files(findex(6)).name); end;
                    else % dicom
                        for ff=1:6
                            if findex(ff)
                                fnames{count}{ff} = fullfile(dataset_folder,ids(i).name,files(findex(ff)).name);
                                dcmfiles = dir(fnames{count}{ff});
                                for dd=1:numel(dcmfiles)
                                    if (~strcmp(dcmfiles(dd).name(1),'.') && ~dcmfiles(dd).isdir)
                                        try
                                            dicominfo(fullfile(fnames{count}{ff},dcmfiles(dd).name));
                                            dflag = 1;
                                        catch
                                            dflag = 0;
                                        end
                                        if dflag
                                            % found dicom file we think
                                            fnames{count}{ff} = fullfile(fnames{count}{ff},dcmfiles(dd).name);
                                            break;
                                        end
                                    end
                                end
                            end
                        end %for
                    end
                end
                
            case 2 % subdirectories contain data
                tpts = dir(fullfile(dataset_folder,ids(i).name));
                
                % cycle through data point subdirectories for an ID
                for j=1:numel(tpts)
                    if (~strcmp(tpts(j).name(1),'.') && tpts(j).isdir)
                        % process if is directory and not a hidden file
                        fprintf('  data point = %s\n',tpts(j).name);
                        cd(fullfile(dataset_folder,ids(i).name,tpts(j).name));
                        %files = dir('*.fld');
                        files = dir();
                        findex = zeros(1,6);
                        for k=1:numel(files)
                            for m=1:6
                                if findf(m) &&...
                                        numel(strfind(files(k).name,f_format{m}))
                                    findex(m) = k;
                                end
                            end
                        end
                        % print out results of search
                        if sum(findex)==0
                            fprintf(' ERROR - 0 files found\n');
                        else
                            % Found some or all
                            for m=1:6
                                if findf(m)
                                    if findex(m)
                                        fprintf('        %s = %s\n',f_format{m},files(findex(m)).name);
                                    else
                                        fprintf('        %s = ERROR no file found.\n',f_format{m});
                                    end
                                end
                            end
                        end
                        
                        if findex(1)>0 && findex(2)>0 && findex(3)>0 % found more than 1 file
                            count = count + 1;
                            if ~dicomFlag
                                fnames{count}{1} = fullfile(dataset_folder,ids(i).name,tpts(j).name,files(findex(1)).name);
                                fnames{count}{2} = fullfile(dataset_folder,ids(i).name,tpts(j).name,files(findex(2)).name);
                                fnames{count}{3} = fullfile(dataset_folder,ids(i).name,tpts(j).name,files(findex(3)).name);
                                if findex(4) fnames{count}{4} = fullfile(dataset_folder,ids(i).name,tpts(j).name,files(findex(4)).name); end;
                                if findex(5) fnames{count}{5} = fullfile(dataset_folder,ids(i).name,tpts(j).name,files(findex(5)).name); end;
                                if findex(6) fnames{count}{6} = fullfile(dataset_folder,ids(i).name,tpts(j).name,files(findex(6)).name); end;
                            else % dicom
                                for ff=1:6
                                    if findex(ff)
                                        fnames{count}{ff} = fullfile(dataset_folder,ids(i).name,tpts(j).name,files(findex(ff)).name);
                                        dcmfiles = dir(fnames{count}{ff});
                                        for dd=1:numel(dcmfiles)
                                            if (~strcmp(dcmfiles(dd).name(1),'.') && ~dcmfiles(dd).isdir)
                                                % found dicom file we think
                                                try
                                                    dicominfo(fullfile(fnames{count}{ff},dcmfiles(dd).name));
                                                    dflag = 1;
                                                catch
                                                    dflag = 0;
                                                end
                                                if dflag
                                                    fnames{count}{ff} = fullfile(fnames{count}{ff},dcmfiles(dd).name);
                                                    break;
                                                end
                                            end
                                        end
                                    end % for ff
                                end %for
                            end
                        end
                    end
                end % for timepoints
            case 3 % subsubdirectories contain data
                tpts = dir(fullfile(dataset_folder,ids(i).name));
                
                % cycle through data point subdirectories for an ID
                for j=1:numel(tpts)
                    if (~strcmp(tpts(j).name(1),'.') && tpts(j).isdir)
                        % process if is directory and not a hidden file
                        fprintf('  data point = %s\n',tpts(j).name);
                        cd(fullfile(dataset_folder,ids(i).name,tpts(j).name));
                        
                        % in timepoint
                        % now go down one more directory before looking for
                        % files
                        dirintpts = dir();
                        for l=1:numel(dirintpts);
                            if (~strcmp(dirintpts(l).name(1),'.') && dirintpts(l).isdir)
                                fprintf('  sub data point = %s\n',dirintpts(l).name);
                                cd(fullfile(dataset_folder,ids(i).name,tpts(j).name,dirintpts(l).name));
                                files = dir();
                                findex = zeros(1,6);
                                for k=1:numel(files)
                                    for m=1:6
                                        if findf(m) &&...
                                                numel(strfind(files(k).name,f_format{m}))
                                            findex(m) = k;
                                        end
                                    end
                                end                              
                                
                                % print out results of search
                                if sum(findex)==0
                                    fprintf(' ERROR - 0 files found\n');
                                else
                                    % Found some or all
                                    for m=1:6
                                        if findf(m)
                                            if findex(m)
                                                fprintf('        %s = %s\n',f_format{m},files(findex(m)).name);
                                            else
                                                fprintf('        %s = ERROR no file found.\n',f_format{m});
                                            end
                                        end
                                    end
                                end
                                
                                if findex(1)>0 && findex(2)>0 && findex(3)>0 % found more than 1 file
                                    count = count + 1;
                                    if ~dicomFlag
                                        fnames{count}{1} = fullfile(dataset_folder,ids(i).name,tpts(j).name,dirintpts(l).name,files(findex(1)).name);
                                        fnames{count}{2} = fullfile(dataset_folder,ids(i).name,tpts(j).name,dirintpts(l).name,files(findex(2)).name);
                                        fnames{count}{3} = fullfile(dataset_folder,ids(i).name,tpts(j).name,dirintpts(l).name,files(findex(3)).name);
                                        if findex(4) fnames{count}{4} = fullfile(dataset_folder,ids(i).name,tpts(j).name,dirintpts(l).name,files(findex(4)).name); end;
                                        if findex(5) fnames{count}{5} = fullfile(dataset_folder,ids(i).name,tpts(j).name,dirintpts(l).name,files(findex(5)).name); end;
                                        if findex(6) fnames{count}{6} = fullfile(dataset_folder,ids(i).name,tpts(j).name,dirintpts(l).name,files(findex(6)).name); end;
                                    else % dicom
                                        for ff=1:6
                                            if findex(ff)
                                                fnames{count}{ff} = fullfile(dataset_folder,ids(i).name,tpts(j).name,files(findex(ff)).name);
                                                dcmfiles = dir(fnames{count}{ff});
                                                for dd=1:numel(dcmfiles)
                                                    if (~strcmp(dcmfiles(dd).name(1),'.') && ~dcmfiles(dd).isdir)
                                                        % found dicom file we think
                                                        try
                                                            dicominfo(fullfile(fnames{count}{ff},dcmfiles(dd).name));
                                                            dflag = 1;
                                                        catch
                                                            dflag = 0;
                                                        end
                                                        if dflag
                                                            fnames{count}{ff} = fullfile(fnames{count}{ff},dcmfiles(dd).name);
                                                            break;
                                                        end
                                                    end
                                                end
                                            end % for ff
                                        end %for
                                    end
                                end
                            end % found file set
                        end % dirs in tmpt
                    end
                end % for timepoints
                
            otherwise
                fprintf('ERROR: no method to process levels %d.\n',levels);
                
        end % switch levels
    end
end


cd(savedir);

% fnames{tt}{1} = fullfile(current_dir, filename_exp);
% fnames{tt}{2} = fullfile(current_dir, filename_insR);
% fnames{tt}{3} = fullfile(current_dir, filename_voi);


end % function

