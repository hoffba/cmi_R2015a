% CMI script
function checkPRMfiles(cmiObj)
%   Input: cmiObj = CMIclass object containing current settings

%% Select study directory to process

dataset_folder = uigetdir('/mnt/cmi/projects/CT Lung/clinical/clinical_longitudinal',...
    'Select directory containing all patient subdirectories:');
if (dataset_folder == 0)
    disp('Error: no study main fld directory specified. Abort')
    return;
end

%% Set constant parsing strings

f1_format = '_Exp';
f2_format = '_Ins_R';
f3_format = '_Exp_VOI';
levels = 1;


prompt = {'First filename format:' 'Second filename format:' 'VOI filename format' 'Directory levels'};
title = 'PRM Filename specs';
defAns = {f1_format f2_format f3_format num2str(levels)};
options.Resize = 'on';
answer = inputdlg(prompt,title,4,defAns,options);

f1_format = answer{1};
f2_format = answer{2};
f3_format = answer{3};
levels = str2double(answer{4});

fprintf('\nParse filenames using f1=%s f2=%s f3=%s down %d levels\n',f1_format,f2_format,f3_format,levels);
fprintf('-------------------------------------------------------------\n');

%% Create subdirectory/file lists from input study summary data

savedir = pwd;
ids = dir(dataset_folder);
count = 0;
pcount = 0;

for i=1:numel(ids)
    % Cycle through IDs in experiment
    if (~strcmp(ids(i).name(1),'.') && ids(i).isdir)
        % process if is directory and not a hidden file
        fprintf('Process ID=%s\n',ids(i).name);
        cd(fullfile(dataset_folder,ids(i).name));
        
        % Process either at this level or in subdirectories
        switch levels
            case 1 % data in this directory
                files = dir('*.fld');
                f1 = 0; f2 = 0; f3 = 0;
                for k=1:numel(files)
                    if numel(strfind(files(k).name,f1_format)) &&...
                            ~numel(strfind(files(k).name,'_VOI'))
                        f1 = k;
                    end
                    if numel(strfind(files(k).name,f2_format)) &&...
                            ~numel(strfind(files(k).name,'VOI'))
                        f2 = k;
                    end
                    if numel(strfind(files(k).name,f3_format)) &&...
                            ~numel(strfind(files(k).name,'VOI_gr'))
                        f3 = k;
                    end
                end
                % print out results of search
                if f1==0 && f2==0 && f3==0
                    fprintf(' ERROR - 0 files found\n');
                else
                    % Found some or all
                    if f1 fprintf('        %s = %s\n',f1_format,files(f1).name);
                    else  fprintf('        %s = ERROR no file found.\n',f1_format); end
                    if f2 fprintf('        %s = %s\n',f2_format,files(f2).name);
                    else  fprintf('        %s = ERROR no file found.\n',f2_format); end
                    if f3 fprintf('        %s = %s\n',f3_format,files(f3).name);
                    else  fprintf('        %s = ERROR no file found.\n',f3_format); end
                    pcount = pcount + 1;
                end
                
                if f1~=0 && f2~=0 && f3~=0
                    count = count + 1;
                    fnames{count}{1} = fullfile(dataset_folder,ids(i).name,files(f1).name);
                    fnames{count}{2} = fullfile(dataset_folder,ids(i).name,files(f2).name);
                    fnames{count}{3} = fullfile(dataset_folder,ids(i).name,files(f3).name);
                end
                
            case 2 % subdirectories contain data
                
                
                tpts = dir(fullfile(dataset_folder,ids(i).name));
                
                % cycle through data points for an ID
                for j=1:numel(tpts)
                    if (~strcmp(tpts(j).name(1),'.') && tpts(j).isdir)
                        % process if is directory and not a hidden file
                        fprintf('  Data point=%s',tpts(j).name);
                        cd(fullfile(fullfile(dataset_folder,ids(i).name,tpts(j).name)));
                        files = dir('*.fld');
                        
                        f1 = 0; f2 = 0; f3 = 0;
                        for k=1:numel(files)
                            if numel(strfind(files(k).name,f1_format)) &&...
                                    ~numel(strfind(files(k).name,'_VOI'))
                                %fprintf('       AVS _Ins field file = %s\n',files(k).name)
                                f1 = k;
                            end
                            if numel(strfind(files(k).name,f2_format)) &&...
                                    ~numel(strfind(files(k).name,'VOI'))
                                %fprintf('       AVS _Exp field file = %s\n',files(k).name)
                                f2 = k;
                            end
                            if numel(strfind(files(k).name,f3_format))
                                %fprintf('       AVS _VOI field file = %s\n',files(k).name)
                                f3 = k;
                            end
                        end
                        
                        % print out results of search
                        if f1~=0 && f2~=0 && f3~=0
                            fprintf(' GOOD - found all 3 files\n');
                        elseif f1==0 && f2==0 && f3==0
                            fprintf(' ERROR - 0 files found\n');
                        else
                            % Found some
                            fprintf('\n');
                            if f1 fprintf('        %s = %s\n',f1_format,files(f1).name);
                            else  fprintf('        %s = ERROR no file found.\n',f1_format); end
                            if f2 fprintf('        %s = %s\n',f2_format,files(f2).name);
                            else  fprintf('        %s = ERROR no file found.\n',f2_format); end
                            if f3 fprintf('        %s = %s\n',f3_format,files(f3).name);
                            else  fprintf('        %s = ERROR no file found.\n',f3_format); end
                            pcount = pcount + 1;
                        end
                        
                        if f1~=0 && f2~=0 && f3~=0
                            count = count + 1;
                            fnames{count}{1} = fullfile(dataset_folder,ids(i).name,tpts(j).name,files(f1).name);
                            fnames{count}{2} = fullfile(dataset_folder,ids(i).name,tpts(j).name,files(f2).name);
                            fnames{count}{3} = fullfile(dataset_folder,ids(i).name,tpts(j).name,files(f3).name);
                        end
                    end
                end
            otherwise
                fprintf('ERROR: no method to process %d levels\n',levels);
        end
    end
end

fprintf('DONE checking PRM files with %d full triplets, %d partial sets.\n',count,pcount);

cd(savedir);