function updateYACTAini

[~,upath] = system('echo %USERPROFILE%');
upath = strtrim(upath);

ypath = fullfile(upath,'AppData','Roaming','yacta64');

%% Update yacta.ini
fn = fullfile(ypath,'yacta.ini');
fprintf('Updating file: %s\n',fn);
str = readlines(fn);
ind = find(startsWith(str,'ExportJobsMHD'),1);
str(ind) = "ExportJobsMHD = lung_lobes.labelexport lung_right_left.labelexport tbt_lobes.labelexport";
ind = find(startsWith(str,'SaveLabelFormat'),1);
str(ind) = "SaveLabelFormat=1";
writelines(str,fn);

%% Write supporting import/export files:
fn = fullfile(ypath,'lung.labelexport');
if ~isfile(fn)
    fprintf('Writing file: %s\n',fn);
    writefile(fn,{{'fileformat','mhd'},...
                  {'0','defaultvalue'},...
                  {'10','RightLung'},...
                  {'20','LeftLung'}});
end

fn = fullfile(ypath,'Lung.labelimport');
if ~isfile(fn)
    fprintf('Writing file: %s\n',fn);
    writefile(fn,{{'fileformat','mhd'},...
                  {'maskdir	D:\MedicalImages\hd_segmentations'},...
                  {'delete','LeftLung','RightLung','Lung'},...
                  {'add','.*LeftLung\.mhd','LeftLung'},...
                  {'add','.*LeftLung\.mhd','Lung'},...
                  {'add','.*RightLung\.mhd','RightLung'},...
                  {'add','.*RightLung\.mhd','Lung'}});
end

fn = fullfile(ypath,'lung_lobes.labelexport');
if ~isfile(fn)
    fprintf('Writing file: %s\n',fn);
    writefile(fn,{{'fileformat','mhd'},...
                  {'0','defaultvalue'},...
                  {'10','UpperLobusRightLung'},...
                  {'20','MidLobusRightLung'},...
                  {'30','LowerLobusRightLung'},...
                  {'40','UpperLobusLeftLung'},...
                  {'50','LingulaLeftLung'},...
                  {'60','LowerLobusLeftLung'}});
end

fn = fullfile(ypath,'lung_right_left.labelexport');
if ~isfile(fn)
    fprintf('Writing file: %s\n',fn);
    writefile(fn,{{'fileformat','mhd'},...
                  {'0','defaultvalue'},...
                  {'10','RightLung'},...
                  {'20','LeftLung'}});
end



function writefile(fn,C)
Ni = numel(C);
fid = fopen(fn,'w');
for i = 1:Ni
    Nj = numel(C{i});
    for j = 1:Nj
        fprintf(fid,C{i}{j});
        if j~=Nj
            fprintf(fid,'\t');
        end
    end
    if i~=Ni
        fprintf(fid,'\n');
    end
end
fclose(fid);