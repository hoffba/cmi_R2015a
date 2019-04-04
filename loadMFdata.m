function [stat,datestrs] = loadMFdata(cmiObj,subjdir)

stat = true;
dnames = dir(subjdir);
dnames = {dnames([dnames.isdir]).name}';
dnames = dnames(cellfun(@length,dnames)==8);
datenums = cellfun(@(x)datenum(x,'yyyymmdd'),dnames);
datenums(datenums==737061) = [];
datestrs = cellfun(@(x)datestr(x,'yyyymmdd'),num2cell(datenums),'UniformOutput',false);

[~,subjID] = fileparts(subjdir);

% Load baseline images and VOI:
tags = {'MT_MTR';'FWI_FatPct';'FWI_R2star';'DWI_T2w';'DWI_HighB';'DWI_ADC'};
rnames = fullfile(subjdir,datestrs{1},...
    strcat(subjID,'_',datestrs{1},'_',[{'MT_off'};tags],'_R.mhd'));
fnames = fullfile(subjdir,datestrs{1},...
    strcat(subjID,'_',datestrs{1},'_',[{'MT_off'};tags],'.mhd'));
ri = logical(cellfun(@(x)exist(x,'file'),rnames));
fnames(ri) = rnames(ri);
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
            datestrs{i} = '';
        end
    end
    
else
    fprintf('\nFiles not found:\n');
    fprintf('%s\n',fnames{~vi});
    stat = false;
    datestrs = {};
end
datestrs = datestrs(~cellfun(@isempty,datestrs));