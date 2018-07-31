function loadMFdata(cmiObj,subjdir)

dnames = dir(subjdir);
dnames = {dnames([dnames.isdir]).name}';
datenums = cellfun(@(x)datenum(x,'yyyymmdd'),dnames);
datenums(datenums==737061) = [];
datestrs = cellfun(@(x)datestr(x,'yyyymmdd'),num2cell(datenums),'UniformOutput',false);

[~,subjID] = fileparts(subjdir);

% Load baseline images and VOI:
tags = {'MT_MTR';'FWI_FatPct';'FWI_R2star';'DWI_T2w';'DWI_HighB';'DWI_ADC'};
fnames = fullfile(subjdir,datestrs{1},...
    strcat(subjID,'_',datestrs{1},'_',[{'MT_off'};tags],'.mhd'));
vi = logical(cellfun(@(x)exist(x,'file'),fnames));
if all(vi)
    
    cmiObj.loadImg(0,fnames);
    
    vname = fullfile(subjdir,datestrs{1},...
        strcat(subjID,'_',datestrs{1},'_MT_VOI.mhd'));
    if exist(vname,'file')
        cmiObj.loadMask(vname);
    else
        fprintf('\nFile not found:\n');
        fprintf('%s\n',vname);
    end
    
    for i = 2:length(datestrs)
        fnames = fullfile(subjdir,datestrs{i},...
            strcat(subjID,'_',datestrs{i},'_',tags,'_R.mhd'));
        vi = logical(cellfun(@(x)exist(x,'file'),fnames));
        if all(vi)
            cmiObj.loadImg(1,fnames);
        else
            fprintf('\nFiles not found:\n');
            fprintf('%s\n',fnames{~vi});
        end
    end
    
else
    fprintf('\nFiles not found:\n');
    fprintf('%s\n',fnames{~vi});
end